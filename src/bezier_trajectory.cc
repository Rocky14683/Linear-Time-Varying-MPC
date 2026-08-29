#include "differential_drive_mpc/bezier_trajectory.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace differential_drive_mpc {
namespace {

constexpr int kConstraintSamples = 10001;
constexpr int kArcLengthSamples = 20001;
constexpr double kConstraintSafetyFactor = 1.000001;
constexpr double kMinimumCurveDerivative = 1e-5;

struct ProgressProfile {
  double value;
  double first_derivative;
  double second_derivative;
};

// A minimum-jerk normalized progress profile. Its first and second derivatives
// vanish at both endpoints, producing a trajectory that starts and ends at rest.
ProgressProfile MinimumJerkProgress(double normalized_time) {
  const double t = std::clamp(normalized_time, 0.0, 1.0);
  const double one_minus_t = 1.0 - t;
  return {.value = t * t * t * (10.0 - 15.0 * t + 6.0 * t * t),
          .first_derivative = 30.0 * t * t * one_minus_t * one_minus_t,
          .second_derivative = 60.0 * t * one_minus_t * (1.0 - 2.0 * t)};
}

Eigen::Vector2d DeCasteljau(std::vector<Eigen::Vector2d> points, double t) {
  for (int degree = static_cast<int>(points.size()) - 1; degree > 0; --degree) {
    for (int index = 0; index < degree; ++index) {
      points[index] = (1.0 - t) * points[index] + t * points[index + 1];
    }
  }
  return points.front();
}

std::vector<Eigen::Vector2d> DerivativeControlPoints(std::vector<Eigen::Vector2d> points,
                                                     int order) {
  for (int derivative = 0; derivative < order; ++derivative) {
    if (points.size() < 2) {
      return {Eigen::Vector2d::Zero()};
    }
    const double degree = static_cast<double>(points.size() - 1);
    std::vector<Eigen::Vector2d> next;
    next.reserve(points.size() - 1);
    for (size_t index = 0; index + 1 < points.size(); ++index) {
      next.push_back(degree * (points[index + 1] - points[index]));
    }
    points = std::move(next);
  }
  return points;
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

BezierTrajectory BezierTrajectory::Create(std::vector<Eigen::Vector2d> control_points,
                                          double requested_duration,
                                          const TrajectoryLimits& limits) {
  if (control_points.size() < 4) {
    throw std::invalid_argument("A cubic-or-higher Bezier curve requires at least four points.");
  }
  if (requested_duration <= 0.0) {
    throw std::invalid_argument("The requested duration must be positive.");
  }
  ValidateLimits(limits);

  BezierTrajectory trajectory(std::move(control_points), requested_duration, limits);
  trajectory.duration_ *= trajectory.FindRequiredTimeScale() * kConstraintSafetyFactor;
  return trajectory;
}

BezierTrajectory::BezierTrajectory(std::vector<Eigen::Vector2d> control_points, double duration,
                                   TrajectoryLimits limits)
    : control_points_(std::move(control_points)), duration_(duration), limits_(limits) {
  BuildArcLengthTable();
}

TrajectoryPoint BezierTrajectory::Evaluate(double time) const {
  const double clamped_time = std::clamp(time, 0.0, duration_);
  const double normalized_time = clamped_time / duration_;
  const ProgressProfile progress = MinimumJerkProgress(normalized_time);
  const double distance = total_length_ * progress.value;
  const double bezier_parameter = ParameterAtDistance(distance);
  const Eigen::Vector2d position = EvaluateDerivative(bezier_parameter, 0);
  const Eigen::Vector2d path_derivative = EvaluateDerivative(bezier_parameter, 1);
  const Eigen::Vector2d path_second_derivative = EvaluateDerivative(bezier_parameter, 2);
  const Eigen::Vector2d path_third_derivative = EvaluateDerivative(bezier_parameter, 3);
  const double path_speed_squared = path_derivative.squaredNorm();
  if (path_speed_squared < kMinimumCurveDerivative * kMinimumCurveDerivative) {
    throw std::runtime_error("Bezier trajectory has a cusp or near-zero path derivative.");
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
  const double speed = total_length_ * progress.first_derivative / duration_;
  const double linear_acceleration =
      total_length_ * progress.second_derivative / (duration_ * duration_);
  const Eigen::Vector2d velocity = tangent * speed;
  const Eigen::Vector2d acceleration =
      tangent * linear_acceleration + normal * curvature * speed * speed;
  const double angular_velocity = curvature * speed;
  const double angular_acceleration =
      curvature_per_distance * speed * speed + curvature * linear_acceleration;

  return {.time = clamped_time,
          .pose = {.position = position,
                   .heading = std::atan2(path_derivative.y(), path_derivative.x())},
          .velocity = velocity,
          .acceleration = acceleration,
          .twist = {.linear = speed, .angular = angular_velocity},
          .acceleration_twist = {.linear = linear_acceleration, .angular = angular_acceleration}};
}

std::vector<TrajectoryPoint> BezierTrajectory::Sample(int sample_count) const {
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

std::vector<Eigen::Vector2d> BezierTrajectory::SamplePositionsByDistance(int sample_count) const {
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

Eigen::Vector2d BezierTrajectory::EvaluateDerivative(double bezier_parameter, int order) const {
  if (order == 0) {
    return DeCasteljau(control_points_, bezier_parameter);
  }
  return DeCasteljau(DerivativeControlPoints(control_points_, order), bezier_parameter);
}

void BezierTrajectory::BuildArcLengthTable() {
  arc_length_table_.clear();
  arc_length_table_.reserve(kArcLengthSamples);

  Eigen::Vector2d previous_position = EvaluateDerivative(0.0, 0);
  total_length_ = 0.0;
  for (int index = 0; index < kArcLengthSamples; ++index) {
    const double bezier_parameter =
        static_cast<double>(index) / static_cast<double>(kArcLengthSamples - 1);
    const Eigen::Vector2d path_derivative = EvaluateDerivative(bezier_parameter, 1);
    if (path_derivative.squaredNorm() < kMinimumCurveDerivative * kMinimumCurveDerivative) {
      throw std::invalid_argument("Control points create a cusp or near-zero path derivative.");
    }

    const Eigen::Vector2d position = EvaluateDerivative(bezier_parameter, 0);
    if (index > 0) {
      total_length_ += (position - previous_position).norm();
    }
    arc_length_table_.push_back({.bezier_parameter = bezier_parameter, .distance = total_length_});
    previous_position = position;
  }

  if (total_length_ < kMinimumCurveDerivative) {
    throw std::invalid_argument("Control points must define a non-zero-length path.");
  }
}

double BezierTrajectory::ParameterAtDistance(double distance) const {
  const double clamped_distance = std::clamp(distance, 0.0, total_length_);
  const auto upper =
      std::lower_bound(arc_length_table_.begin(), arc_length_table_.end(), clamped_distance,
                       [](const ArcLengthSample& sample, double requested_distance) {
                         return sample.distance < requested_distance;
                       });
  if (upper == arc_length_table_.begin()) {
    return upper->bezier_parameter;
  }
  if (upper == arc_length_table_.end()) {
    return arc_length_table_.back().bezier_parameter;
  }

  const ArcLengthSample& lower_sample = *std::prev(upper);
  const double segment_length = upper->distance - lower_sample.distance;
  if (segment_length <= 0.0) {
    return upper->bezier_parameter;
  }
  const double fraction = (clamped_distance - lower_sample.distance) / segment_length;
  return lower_sample.bezier_parameter +
         fraction * (upper->bezier_parameter - lower_sample.bezier_parameter);
}

double BezierTrajectory::FindRequiredTimeScale() const {
  double required_scale = 1.0;
  for (int index = 0; index < kConstraintSamples; ++index) {
    const double normalized_time = static_cast<double>(index) / (kConstraintSamples - 1);
    const ProgressProfile progress = MinimumJerkProgress(normalized_time);
    const double distance = total_length_ * progress.value;
    const double bezier_parameter = ParameterAtDistance(distance);
    const Eigen::Vector2d path_derivative = EvaluateDerivative(bezier_parameter, 1);
    const Eigen::Vector2d path_second_derivative = EvaluateDerivative(bezier_parameter, 2);
    const Eigen::Vector2d path_third_derivative = EvaluateDerivative(bezier_parameter, 3);
    const double path_speed_squared = path_derivative.squaredNorm();
    if (path_speed_squared < kMinimumCurveDerivative * kMinimumCurveDerivative) {
      throw std::invalid_argument("Control points create a cusp or near-zero path derivative.");
    }
    const double path_speed = std::sqrt(path_speed_squared);
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
    const double speed = total_length_ * progress.first_derivative / duration_;
    const double linear_acceleration =
        total_length_ * progress.second_derivative / (duration_ * duration_);
    const double angular_velocity = curvature * speed;
    const double angular_acceleration =
        curvature_per_distance * speed * speed + curvature * linear_acceleration;
    const double half_track_width = limits_.track_width / 2.0;
    const double left_wheel_velocity =
        (speed - angular_velocity * half_track_width) / limits_.wheel_radius;
    const double right_wheel_velocity =
        (speed + angular_velocity * half_track_width) / limits_.wheel_radius;
    const double left_wheel_acceleration =
        (linear_acceleration - angular_acceleration * half_track_width) / limits_.wheel_radius;
    const double right_wheel_acceleration =
        (linear_acceleration + angular_acceleration * half_track_width) / limits_.wheel_radius;

    required_scale = std::max(required_scale, speed / limits_.max_linear_velocity);
    required_scale = std::max(
        required_scale, std::sqrt(std::abs(linear_acceleration) / limits_.max_linear_acceleration));
    required_scale =
        std::max(required_scale, std::abs(angular_velocity) / limits_.max_angular_velocity);
    required_scale = std::max(required_scale, std::sqrt(std::abs(angular_acceleration) /
                                                        limits_.max_angular_acceleration));
    required_scale =
        std::max(required_scale, std::abs(left_wheel_velocity) / limits_.max_wheel_velocity);
    required_scale =
        std::max(required_scale, std::abs(right_wheel_velocity) / limits_.max_wheel_velocity);
    required_scale = std::max(required_scale, std::sqrt(std::abs(left_wheel_acceleration) /
                                                        limits_.max_wheel_acceleration));
    required_scale = std::max(required_scale, std::sqrt(std::abs(right_wheel_acceleration) /
                                                        limits_.max_wheel_acceleration));
  }
  return required_scale;
}

}  // namespace differential_drive_mpc
