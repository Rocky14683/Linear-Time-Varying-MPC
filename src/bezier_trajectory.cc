#include "differential_drive_mpc/bezier_trajectory.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace differential_drive_mpc {
namespace {

constexpr int kArcLengthSamples = 20001;
// Match the dense verifier's spatial resolution. This is especially important
// when an endpoint speed is prescribed, because a final uniform time stretch
// would no longer preserve that boundary condition.
constexpr int kTimeProfileSamples = 10001;
constexpr double kMinimumCurveDerivative = 1e-5;
constexpr double kConstraintTolerance = 1e-8;
constexpr int kProfileBinarySearchIterations = 50;
constexpr int kProfileVerificationSamples = 10001;

double BasisFunction(int index, int degree, double parameter, const std::vector<double>& knots,
                     int last_control_index) {
  if (degree == 0) {
    return ((knots[index] <= parameter && parameter < knots[index + 1]) ||
            (parameter == knots.back() && index == last_control_index))
               ? 1.0
               : 0.0;
  }
  double value = 0.0;
  const double left_denominator = knots[index + degree] - knots[index];
  if (left_denominator > 0.0) {
    value += (parameter - knots[index]) / left_denominator *
             BasisFunction(index, degree - 1, parameter, knots, last_control_index);
  }
  const double right_denominator = knots[index + degree + 1] - knots[index + 1];
  if (right_denominator > 0.0) {
    value += (knots[index + degree + 1] - parameter) / right_denominator *
             BasisFunction(index + 1, degree - 1, parameter, knots, last_control_index);
  }
  return value;
}

Eigen::Vector2d EvaluateBSpline(const std::vector<Eigen::Vector2d>& control_points,
                                const std::vector<double>& knots, int degree, double parameter) {
  const int last_control_index = static_cast<int>(control_points.size()) - 1;
  const double clamped_parameter =
      std::clamp(parameter, knots[degree], knots[last_control_index + 1]);
  Eigen::Vector2d result = Eigen::Vector2d::Zero();
  for (int index = 0; index <= last_control_index; ++index) {
    result += BasisFunction(index, degree, clamped_parameter, knots, last_control_index) *
              control_points[index];
  }
  return result;
}

void ValidateLimits(const TrajectoryLimits& limits) {
  if (limits.max_linear_velocity <= 0.0 || limits.max_linear_acceleration <= 0.0 ||
      limits.max_angular_velocity <= 0.0 || limits.max_angular_acceleration <= 0.0 ||
      limits.wheel_radius <= 0.0 || limits.track_width <= 0.0 || limits.max_wheel_velocity <= 0.0 ||
      limits.max_wheel_acceleration <= 0.0) {
    throw std::invalid_argument("All trajectory limits must be positive.");
  }
}

}  // namespace

CubicBSplineTrajectory CubicBSplineTrajectory::Create(std::vector<Eigen::Vector2d> waypoints,
                                                       const TrajectoryLimits& limits, bool reverse) {
  if (waypoints.size() < 4) {
    throw std::invalid_argument("A cubic interpolating B-spline requires at least four waypoints.");
  }
  ValidateLimits(limits);

  CubicBSplineTrajectory trajectory(std::move(waypoints), limits, reverse);
  trajectory.BuildMinimumTimeProfile();
  return trajectory;
}

CubicBSplineTrajectory::CubicBSplineTrajectory(std::vector<Eigen::Vector2d> waypoints,
                                                TrajectoryLimits limits, bool reverse)
    : waypoints_(std::move(waypoints)),
      limits_(limits),
      reverse_(reverse) {
  BuildInterpolatingSpline();
  BuildArcLengthTable();
}

TrajectoryPoint CubicBSplineTrajectory::EvaluateImpl(double time) const {
  const double clamped_time = std::clamp(time, 0.0, duration_);
  const ProgressState progress = ProgressAtTime(clamped_time);
  const double distance = progress.distance;
  const double path_parameter = ParameterAtDistance(distance);
  const Eigen::Vector2d position = EvaluateDerivative(path_parameter, 0);
  const Eigen::Vector2d path_derivative = EvaluateDerivative(path_parameter, 1);
  const Eigen::Vector2d path_second_derivative = EvaluateDerivative(path_parameter, 2);
  const Eigen::Vector2d path_third_derivative = EvaluateDerivative(path_parameter, 3);
  const double path_speed_squared = path_derivative.squaredNorm();
  if (path_speed_squared < kMinimumCurveDerivative * kMinimumCurveDerivative) {
    throw std::runtime_error("B-spline trajectory has a cusp or near-zero path derivative.");
  }

  const double path_speed = std::sqrt(path_speed_squared);
  const Eigen::Vector2d tangent = path_derivative / path_speed;
  const Eigen::Vector2d normal(-tangent.y(), tangent.x());
  const double path_cross = path_derivative.x() * path_second_derivative.y() -
                            path_derivative.y() * path_second_derivative.x();
  const double curvature = path_cross / (path_speed * path_speed_squared);
  const double path_cross_derivative = path_derivative.x() * path_third_derivative.y() -
                                       path_derivative.y() * path_third_derivative.x();
  const double path_speed_derivative = path_derivative.dot(path_second_derivative) / path_speed;
  const double curvature_per_distance =
      path_cross_derivative / (path_speed_squared * path_speed_squared) -
      3.0 * path_cross * path_speed_derivative /
          (path_speed_squared * path_speed_squared * path_speed);
  const double speed = progress.speed;
  const double linear_acceleration = progress.acceleration;
  const Eigen::Vector2d velocity = tangent * speed;
  const Eigen::Vector2d acceleration =
      tangent * linear_acceleration + normal * curvature * speed * speed;
  const double angular_velocity = curvature * speed;
  const double angular_acceleration =
      curvature_per_distance * speed * speed + curvature * linear_acceleration;
  const double tangent_heading = std::atan2(path_derivative.y(), path_derivative.x());
  const double heading = reverse_ ? WrapAngle(tangent_heading + kPi) : tangent_heading;
  const double signed_linear_speed = reverse_ ? -speed : speed;
  const double signed_linear_acceleration = reverse_ ? -linear_acceleration : linear_acceleration;

  return {.time = clamped_time,
          .pose = {.position = position,
                   .heading = heading},
          .velocity = velocity,
          .acceleration = acceleration,
          .twist = {.linear = signed_linear_speed, .angular = angular_velocity},
          .acceleration_twist = {.linear = signed_linear_acceleration,
                                 .angular = angular_acceleration}};
}

std::vector<TrajectoryPoint> CubicBSplineTrajectory::Sample(int sample_count) const {
  if (sample_count < 2) {
    throw std::invalid_argument("At least two trajectory samples are required.");
  }
  std::vector<TrajectoryPoint> samples;
  samples.reserve(sample_count);
  for (int index = 0; index < sample_count; ++index) {
    const double fraction = static_cast<double>(index) / (sample_count - 1);
    samples.push_back(Evaluate(fraction * duration_));
  }
  return samples;
}

std::vector<Eigen::Vector2d> CubicBSplineTrajectory::SamplePositionsByDistance(int sample_count) const {
  if (sample_count < 2) {
    throw std::invalid_argument("At least two trajectory samples are required.");
  }

  std::vector<Eigen::Vector2d> samples;
  samples.reserve(sample_count);
  for (int index = 0; index < sample_count; ++index) {
    const double fraction = static_cast<double>(index) / (sample_count - 1);
    const double distance = fraction * total_length_;
    samples.push_back(EvaluateDerivative(ParameterAtDistance(distance), 0));
  }
  return samples;
}

void CubicBSplineTrajectory::BuildInterpolatingSpline() {
  constexpr int kDegree = 3;
  const int waypoint_count = static_cast<int>(waypoints_.size());
  std::vector<double> parameters(waypoint_count, 0.0);
  double total_chord_length = 0.0;
  for (int index = 1; index < waypoint_count; ++index) {
    const double segment_length = (waypoints_[index] - waypoints_[index - 1]).norm();
    if (segment_length < kMinimumCurveDerivative) {
      throw std::invalid_argument("Consecutive waypoints must be distinct.");
    }
    total_chord_length += segment_length;
    parameters[index] = total_chord_length;
  }
  for (double& parameter : parameters) {
    parameter /= total_chord_length;
  }

  knots_.assign(waypoint_count + kDegree + 1, 0.0);
  for (int index = waypoint_count; index < waypoint_count + kDegree + 1; ++index) {
    knots_[index] = 1.0;
  }
  for (int index = 1; index < waypoint_count - kDegree; ++index) {
    knots_[index + kDegree] =
        (parameters[index] + parameters[index + 1] + parameters[index + 2]) / kDegree;
  }

  Eigen::MatrixXd interpolation = Eigen::MatrixXd::Zero(waypoint_count, waypoint_count);
  Eigen::MatrixXd positions(waypoint_count, 2);
  for (int row = 0; row < waypoint_count; ++row) {
    positions.row(row) = waypoints_[row].transpose();
    for (int column = 0; column < waypoint_count; ++column) {
      interpolation(row, column) =
          BasisFunction(column, kDegree, parameters[row], knots_, waypoint_count - 1);
    }
  }
  const Eigen::MatrixXd solved_controls = interpolation.colPivHouseholderQr().solve(positions);
  if ((interpolation * solved_controls - positions).norm() > 1e-9) {
    throw std::runtime_error("Failed to construct an interpolating cubic B-spline.");
  }
  spline_control_points_.resize(waypoint_count);
  for (int index = 0; index < waypoint_count; ++index) {
    spline_control_points_[index] = solved_controls.row(index).transpose();
  }
}

Eigen::Vector2d CubicBSplineTrajectory::EvaluateDerivative(double path_parameter, int order) const {
  std::vector<Eigen::Vector2d> derivative_controls = spline_control_points_;
  std::vector<double> derivative_knots = knots_;
  int degree = 3;
  for (int derivative = 0; derivative < order; ++derivative) {
    if (degree == 0 || derivative_controls.size() < 2) {
      return Eigen::Vector2d::Zero();
    }
    std::vector<Eigen::Vector2d> next_controls;
    next_controls.reserve(derivative_controls.size() - 1);
    for (size_t index = 0; index + 1 < derivative_controls.size(); ++index) {
      const double denominator = derivative_knots[index + degree + 1] - derivative_knots[index + 1];
      next_controls.push_back(degree * (derivative_controls[index + 1] - derivative_controls[index]) /
                              denominator);
    }
    derivative_controls = std::move(next_controls);
    derivative_knots.erase(derivative_knots.begin());
    derivative_knots.pop_back();
    --degree;
  }
  return EvaluateBSpline(derivative_controls, derivative_knots, degree, path_parameter);
}

void CubicBSplineTrajectory::BuildArcLengthTable() {
  arc_length_table_.clear();
  arc_length_table_.reserve(kArcLengthSamples);

  Eigen::Vector2d previous_position = EvaluateDerivative(0.0, 0);
  total_length_ = 0.0;
  for (int index = 0; index < kArcLengthSamples; ++index) {
    const double path_parameter =
        static_cast<double>(index) / static_cast<double>(kArcLengthSamples - 1);
    const Eigen::Vector2d path_derivative = EvaluateDerivative(path_parameter, 1);
    if (path_derivative.squaredNorm() < kMinimumCurveDerivative * kMinimumCurveDerivative) {
      throw std::invalid_argument("Control points create a cusp or near-zero path derivative.");
    }

    const Eigen::Vector2d position = EvaluateDerivative(path_parameter, 0);
    if (index > 0) {
      total_length_ += (position - previous_position).norm();
    }
    arc_length_table_.push_back({.path_parameter = path_parameter, .distance = total_length_});
    previous_position = position;
  }

  if (total_length_ < kMinimumCurveDerivative) {
    throw std::invalid_argument("Control points must define a non-zero-length path.");
  }
}

double CubicBSplineTrajectory::ParameterAtDistance(double distance) const {
  const double clamped_distance = std::clamp(distance, 0.0, total_length_);
  const auto upper =
      std::lower_bound(arc_length_table_.begin(), arc_length_table_.end(), clamped_distance,
                       [](const ArcLengthSample& sample, double requested_distance) {
                         return sample.distance < requested_distance;
                       });
  if (upper == arc_length_table_.begin()) {
    return upper->path_parameter;
  }
  if (upper == arc_length_table_.end()) {
    return arc_length_table_.back().path_parameter;
  }

  const ArcLengthSample& lower_sample = *std::prev(upper);
  const double segment_length = upper->distance - lower_sample.distance;
  if (segment_length <= 0.0) {
    return upper->path_parameter;
  }
  const double fraction = (clamped_distance - lower_sample.distance) / segment_length;
  return lower_sample.path_parameter +
         fraction * (upper->path_parameter - lower_sample.path_parameter);
}

double CubicBSplineTrajectory::CurvatureAtDistance(double distance) const {
  const double parameter = ParameterAtDistance(distance);
  const Eigen::Vector2d first = EvaluateDerivative(parameter, 1);
  const Eigen::Vector2d second = EvaluateDerivative(parameter, 2);
  const double speed_squared = first.squaredNorm();
  const double speed = std::sqrt(speed_squared);
  return (first.x() * second.y() - first.y() * second.x()) / (speed * speed_squared);
}

double CubicBSplineTrajectory::CurvatureDerivativeAtDistance(double distance) const {
  const double parameter = ParameterAtDistance(distance);
  const Eigen::Vector2d first = EvaluateDerivative(parameter, 1);
  const Eigen::Vector2d second = EvaluateDerivative(parameter, 2);
  const Eigen::Vector2d third = EvaluateDerivative(parameter, 3);
  const double speed_squared = first.squaredNorm();
  const double speed = std::sqrt(speed_squared);
  const double cross = first.x() * second.y() - first.y() * second.x();
  const double cross_derivative = first.x() * third.y() - first.y() * third.x();
  const double speed_derivative = first.dot(second) / speed;
  return cross_derivative / (speed_squared * speed_squared) -
         3.0 * cross * speed_derivative / (speed_squared * speed_squared * speed);
}

bool CubicBSplineTrajectory::AccelerationBoundsAtDistance(double distance, double speed,
                                                           double* lower, double* upper) const {
  *lower = -limits_.max_linear_acceleration;
  *upper = limits_.max_linear_acceleration;
  const double curvature = CurvatureAtDistance(distance);
  const double curvature_derivative = CurvatureDerivativeAtDistance(distance);
  const double curvature_rate = curvature_derivative * speed * speed;

  const auto add_bounded_affine_constraint = [&](double coefficient, double offset,
                                                  double limit) {
    if (std::abs(coefficient) < 1e-10) {
      return std::abs(offset) <= limit + kConstraintTolerance;
    }
    const double first = (-limit - offset) / coefficient;
    const double second = (limit - offset) / coefficient;
    *lower = std::max(*lower, std::min(first, second));
    *upper = std::min(*upper, std::max(first, second));
    return *lower <= *upper + kConstraintTolerance;
  };

  if (!add_bounded_affine_constraint(curvature, curvature_rate,
                                     limits_.max_angular_acceleration)) {
    return false;
  }

  const double half_track = limits_.track_width / 2.0;
  const double wheel_acceleration_limit = limits_.wheel_radius * limits_.max_wheel_acceleration;
  if (!add_bounded_affine_constraint(1.0 - half_track * curvature,
                                     -half_track * curvature_rate,
                                     wheel_acceleration_limit) ||
      !add_bounded_affine_constraint(1.0 + half_track * curvature,
                                     half_track * curvature_rate,
                                     wheel_acceleration_limit)) {
    return false;
  }
  return true;
}

double CubicBSplineTrajectory::VelocityLimitAtDistance(double distance) const {
  const double curvature = CurvatureAtDistance(distance);
  double limit = limits_.max_linear_velocity;
  if (std::abs(curvature) > 1e-10) {
    limit = std::min(limit, limits_.max_angular_velocity / std::abs(curvature));
  }
  const double half_track = limits_.track_width / 2.0;
  limit = std::min(limit, limits_.wheel_radius * limits_.max_wheel_velocity /
                              std::abs(1.0 - half_track * curvature));
  limit = std::min(limit, limits_.wheel_radius * limits_.max_wheel_velocity /
                              std::abs(1.0 + half_track * curvature));

  // Curvature-rate terms consume angular and wheel acceleration even when the
  // path-speed acceleration is zero. Find the largest speed with a non-empty
  // feasible acceleration interval.
  double lower = 0.0;
  double upper = limit;
  for (int iteration = 0; iteration < kProfileBinarySearchIterations; ++iteration) {
    const double candidate = (lower + upper) / 2.0;
    double acceleration_lower = 0.0;
    double acceleration_upper = 0.0;
    if (AccelerationBoundsAtDistance(distance, candidate, &acceleration_lower,
                                     &acceleration_upper)) {
      lower = candidate;
    } else {
      upper = candidate;
    }
  }
  return lower;
}

void CubicBSplineTrajectory::BuildMinimumTimeProfile() {
  const double spacing = total_length_ / (kTimeProfileSamples - 1);
  std::vector<double> velocity_limits(kTimeProfileSamples);
  std::vector<double> maximum_accelerations(kTimeProfileSamples);
  std::vector<double> maximum_decelerations(kTimeProfileSamples);
  for (int index = 0; index < kTimeProfileSamples; ++index) {
    const double distance = index * spacing;
    // The profile is evaluated between spatial nodes, so reserve a small
    // margin for interpolation and the changing curvature within a segment.
    const double physical_velocity_limit = VelocityLimitAtDistance(distance);
    velocity_limits[index] = 0.99 * physical_velocity_limit;
    double lower = 0.0;
    double upper = 0.0;
    // Curvature-rate terms can leave only one signed path acceleration
    // feasible near the velocity cap. Lower the envelope speed until both
    // acceleration and braking remain available; a differential drive can
    // always recover feasibility by slowing down on a regular path.
    double envelope_speed = velocity_limits[index];
    bool has_two_sided_acceleration = false;
    for (int attempt = 0; attempt < 24; ++attempt) {
      if (AccelerationBoundsAtDistance(distance, envelope_speed, &lower, &upper) &&
          lower <= 0.0 && upper >= 0.0) {
        has_two_sided_acceleration = true;
        break;
      }
      envelope_speed *= 0.5;
    }
    if (!has_two_sided_acceleration) {
      throw std::runtime_error("No dynamically feasible time scaling exists for this path.");
    }
    // The two-sided acceleration envelope is also a velocity constraint. If
    // we kept the larger curvature-only cap here, a later pass could speed
    // back into a region where neither braking nor acceleration is feasible.
    if (index != 0 && index != kTimeProfileSamples - 1) {
      velocity_limits[index] = envelope_speed;
    }
    maximum_accelerations[index] = std::max(0.0, upper);
    maximum_decelerations[index] = std::max(0.0, -lower);
  }

  // Forward envelope from rest.
  std::vector<double> forward_speeds(kTimeProfileSamples, 0.0);
  for (int index = 0; index + 1 < kTimeProfileSamples; ++index) {
    const double reachable_speed =
        std::sqrt(forward_speeds[index] * forward_speeds[index] +
                  2.0 * maximum_accelerations[index] * spacing);
    forward_speeds[index + 1] = std::min(velocity_limits[index + 1], reachable_speed);
  }

  // Backward envelope to rest.
  std::vector<double> backward_speeds(kTimeProfileSamples, 0.0);
  for (int index = kTimeProfileSamples - 2; index >= 0; --index) {
    const double brake_limited_speed = std::sqrt(backward_speeds[index + 1] *
                                                     backward_speeds[index + 1] +
                                                 2.0 * maximum_decelerations[index] * spacing);
    backward_speeds[index] = std::min(velocity_limits[index], brake_limited_speed);
  }
  std::vector<double> speeds(kTimeProfileSamples, 0.0);
  for (int index = 0; index < kTimeProfileSamples; ++index) {
    speeds[index] = std::min(forward_speeds[index], backward_speeds[index]);
  }

  time_profile_.clear();
  time_profile_.reserve(kTimeProfileSamples);
  time_profile_.push_back({.distance = 0.0, .time = 0.0, .speed = 0.0});
  for (int index = 0; index + 1 < kTimeProfileSamples; ++index) {
    const double average_speed = (speeds[index] + speeds[index + 1]) / 2.0;
    if (average_speed <= 1e-10) {
      throw std::runtime_error("No dynamically feasible time scaling exists for this path.");
    }
    const TimeProfileSample& previous = time_profile_.back();
    time_profile_.push_back({.distance = (index + 1) * spacing,
                             .time = previous.time + spacing / average_speed,
                             .speed = speeds[index + 1]});
  }

  duration_ = time_profile_.back().time;

  // The profile is generated on a finite spatial grid. Verify all composed
  // chassis and wheel quantities more densely, then apply only the small
  // global stretch required to cover between-node curvature extrema.
  double verification_scale = 1.0;
  for (int index = 0; index < kProfileVerificationSamples; ++index) {
    const double time = duration_ * index / (kProfileVerificationSamples - 1);
    const ProgressState progress = ProgressAtTime(time);
    const double curvature = CurvatureAtDistance(progress.distance);
    const double curvature_derivative = CurvatureDerivativeAtDistance(progress.distance);
    const double angular_velocity = curvature * progress.speed;
    const double angular_acceleration = curvature * progress.acceleration +
                                        curvature_derivative * progress.speed * progress.speed;
    const double half_track = limits_.track_width / 2.0;
    const double left_wheel_velocity =
        (progress.speed - half_track * angular_velocity) / limits_.wheel_radius;
    const double right_wheel_velocity =
        (progress.speed + half_track * angular_velocity) / limits_.wheel_radius;
    const double left_wheel_acceleration =
        (progress.acceleration - half_track * angular_acceleration) / limits_.wheel_radius;
    const double right_wheel_acceleration =
        (progress.acceleration + half_track * angular_acceleration) / limits_.wheel_radius;

    verification_scale = std::max(verification_scale,
                                  progress.speed / limits_.max_linear_velocity);
    verification_scale = std::max(
        verification_scale,
        std::sqrt(std::abs(progress.acceleration) / limits_.max_linear_acceleration));
    verification_scale = std::max(verification_scale,
                                  std::abs(angular_velocity) / limits_.max_angular_velocity);
    verification_scale = std::max(
        verification_scale,
        std::sqrt(std::abs(angular_acceleration) / limits_.max_angular_acceleration));
    verification_scale = std::max(verification_scale,
                                  std::abs(left_wheel_velocity) / limits_.max_wheel_velocity);
    verification_scale = std::max(verification_scale,
                                  std::abs(right_wheel_velocity) / limits_.max_wheel_velocity);
    verification_scale = std::max(
        verification_scale,
        std::sqrt(std::abs(left_wheel_acceleration) / limits_.max_wheel_acceleration));
    verification_scale = std::max(
        verification_scale,
        std::sqrt(std::abs(right_wheel_acceleration) / limits_.max_wheel_acceleration));
  }
  if (verification_scale > 1.0) {
    verification_scale *= 1.001;
    for (TimeProfileSample& sample : time_profile_) {
      sample.time *= verification_scale;
      sample.speed /= verification_scale;
    }
    duration_ = time_profile_.back().time;
  }
}

CubicBSplineTrajectory::ProgressState CubicBSplineTrajectory::ProgressAtTime(double time) const {
  if (time <= 0.0) {
    return {.distance = 0.0, .speed = 0.0, .acceleration = 0.0};
  }
  if (time >= duration_) {
    return {.distance = total_length_, .speed = 0.0, .acceleration = 0.0};
  }
  const auto upper = std::lower_bound(time_profile_.begin(), time_profile_.end(), time,
                                      [](const TimeProfileSample& sample, double requested_time) {
                                        return sample.time < requested_time;
                                      });
  if (upper == time_profile_.begin()) {
    return {.distance = 0.0, .speed = 0.0, .acceleration = 0.0};
  }
  if (upper == time_profile_.end()) {
    return {.distance = total_length_, .speed = 0.0, .acceleration = 0.0};
  }
  const TimeProfileSample& upper_sample = *upper;
  const TimeProfileSample& lower_sample = *std::prev(upper);
  const double interval = upper_sample.time - lower_sample.time;
  const double elapsed = time - lower_sample.time;
  const double acceleration = (upper_sample.speed - lower_sample.speed) / interval;
  return {.distance = std::clamp(lower_sample.distance + lower_sample.speed * elapsed +
                                     0.5 * acceleration * elapsed * elapsed,
                                 lower_sample.distance, upper_sample.distance),
          .speed = lower_sample.speed + acceleration * elapsed,
          .acceleration = acceleration};
}

}  // namespace differential_drive_mpc
