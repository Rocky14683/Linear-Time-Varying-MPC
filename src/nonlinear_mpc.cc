#include "differential_drive_mpc/nonlinear_mpc.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace differential_drive_mpc {
namespace {

Eigen::Vector2d ToVector(Twist2d twist) { return {twist.linear, twist.angular}; }

}  // namespace

LtvMpc::LtvMpc(MpcConfig config) : config_(config) {
  if (config_.horizon_steps < 1 || config_.dt <= 0.0 || config_.optimization_iterations < 1 ||
      config_.gradient_step_size <= 0.0 || config_.max_linear_velocity <= 0.0 ||
      config_.max_angular_velocity <= 0.0 || config_.max_linear_acceleration <= 0.0 ||
      config_.max_angular_acceleration <= 0.0 || config_.wheel_radius <= 0.0 ||
      config_.track_width <= 0.0 || config_.max_wheel_velocity <= 0.0 ||
      config_.max_wheel_acceleration <= 0.0 || config_.linear_response_rate < 0.0 ||
      config_.angular_response_rate < 0.0 || config_.lateral_drift_coefficient < 0.0 ||
      config_.lateral_drift_coefficient > 1.0 || config_.lateral_slip_response_rate < 0.0 ||
      config_.max_traction_acceleration <= 0.0 || config_.linearization_epsilon <= 0.0) {
    throw std::invalid_argument("Invalid dynamic-MPC timing, constraint, or model configuration.");
  }
}

Twist2d LtvMpc::Solve(const RobotState& state, const BezierTrajectory& reference,
                      double reference_time) {
  Eigen::MatrixXd hessian;
  Eigen::VectorXd gradient;
  BuildQuadratic(state, reference, reference_time, &hessian, &gradient);

  ControlSequence controls = InitialGuess(state);
  ProjectControls(&controls, last_command_);
  double best_cost = QuadraticCost(controls, hessian, gradient);

  for (int iteration = 0; iteration < config_.optimization_iterations; ++iteration) {
    Eigen::VectorXd decision(2 * config_.horizon_steps);
    for (int index = 0; index < config_.horizon_steps; ++index) {
      decision.segment<2>(2 * index) = ToVector(controls[index]);
    }
    const Eigen::VectorXd direction = 2.0 * (hessian * decision + gradient);
    if (direction.norm() < 1e-9) {
      break;
    }

    double step_size = config_.gradient_step_size;
    bool accepted = false;
    for (int attempt = 0; attempt < 12; ++attempt) {
      ControlSequence candidate = controls;
      for (int index = 0; index < config_.horizon_steps; ++index) {
        const Eigen::Vector2d update = step_size * direction.segment<2>(2 * index);
        candidate[index].linear -= update.x();
        candidate[index].angular -= update.y();
      }
      ProjectControls(&candidate, last_command_);
      const double candidate_cost = QuadraticCost(candidate, hessian, gradient);
      if (candidate_cost + 1e-10 < best_cost) {
        controls = std::move(candidate);
        best_cost = candidate_cost;
        accepted = true;
        break;
      }
      step_size *= 0.5;
    }
    if (!accepted) {
      break;
    }
  }

  last_command_ = controls.front();
  previous_solution_ = controls;
  predicted_states_ = Rollout(state, controls);
  last_cost_ = best_cost;
  return last_command_;
}

LtvMpc::ControlSequence LtvMpc::InitialGuess(const RobotState& state) {
  ControlSequence controls(config_.horizon_steps, last_command_);
  if (!previous_solution_.empty()) {
    for (int index = 0; index < config_.horizon_steps - 1; ++index) {
      controls[index] = previous_solution_[index + 1];
    }
    controls.back() = previous_solution_.back();
  } else {
    const DynamicState dynamic_state = ToDynamicState(state);
    controls.front() = {.linear = dynamic_state.segment<2>(3).dot(Eigen::Vector2d(
                            std::cos(dynamic_state(2)), std::sin(dynamic_state(2)))),
                        .angular = dynamic_state(5)};
  }
  return controls;
}

void LtvMpc::ProjectControls(ControlSequence* controls, Twist2d previous_control) const {
  const double max_linear_change = config_.max_linear_acceleration * config_.dt;
  const double max_angular_change = config_.max_angular_acceleration * config_.dt;
  const double max_wheel_change = config_.max_wheel_acceleration * config_.dt;
  const double half_track_width = config_.track_width / 2.0;

  for (Twist2d& control : *controls) {
    control.linear =
        std::clamp(control.linear, -config_.max_linear_velocity, config_.max_linear_velocity);
    control.angular =
        std::clamp(control.angular, -config_.max_angular_velocity, config_.max_angular_velocity);
    control.linear = std::clamp(control.linear, previous_control.linear - max_linear_change,
                                previous_control.linear + max_linear_change);
    control.angular = std::clamp(control.angular, previous_control.angular - max_angular_change,
                                 previous_control.angular + max_angular_change);

    double left = (control.linear - control.angular * half_track_width) / config_.wheel_radius;
    double right = (control.linear + control.angular * half_track_width) / config_.wheel_radius;
    const double previous_left =
        (previous_control.linear - previous_control.angular * half_track_width) /
        config_.wheel_radius;
    const double previous_right =
        (previous_control.linear + previous_control.angular * half_track_width) /
        config_.wheel_radius;
    left = std::clamp(left, previous_left - max_wheel_change, previous_left + max_wheel_change);
    right = std::clamp(right, previous_right - max_wheel_change, previous_right + max_wheel_change);
    left = std::clamp(left, -config_.max_wheel_velocity, config_.max_wheel_velocity);
    right = std::clamp(right, -config_.max_wheel_velocity, config_.max_wheel_velocity);
    control.linear = config_.wheel_radius * (left + right) / 2.0;
    control.angular = config_.wheel_radius * (right - left) / config_.track_width;
    previous_control = control;
  }
}

std::vector<RobotState> LtvMpc::Rollout(const RobotState& state,
                                        const ControlSequence& controls) const {
  std::vector<RobotState> states;
  states.reserve(controls.size());
  DynamicState predicted = ToDynamicState(state);
  for (const Twist2d& control : controls) {
    predicted = StepModel(predicted, ToVector(control));
    states.push_back(ToRobotState(predicted));
  }
  return states;
}

LtvMpc::DynamicState LtvMpc::ToDynamicState(const RobotState& state) const {
  Eigen::Vector2d world_velocity = state.world_velocity;
  if (world_velocity.isZero() && std::abs(state.linear_velocity) > 1e-12) {
    world_velocity = Eigen::Vector2d(std::cos(state.pose.heading), std::sin(state.pose.heading)) *
                     state.linear_velocity;
  }
  return {state.pose.position.x(), state.pose.position.y(), state.pose.heading,
          world_velocity.x(),      world_velocity.y(),      state.angular_velocity};
}

RobotState LtvMpc::ToRobotState(const DynamicState& state) const {
  const Eigen::Vector2d heading(std::cos(state(2)), std::sin(state(2)));
  return {.pose = {.position = state.segment<2>(0), .heading = WrapAngle(state(2))},
          .world_velocity = state.segment<2>(3),
          .linear_velocity = state.segment<2>(3).dot(heading),
          .angular_velocity = state(5)};
}

LtvMpc::DynamicState LtvMpc::StepModel(const DynamicState& state,
                                       const Eigen::Vector2d& control) const {
  const double heading_angle = state(2);
  const Eigen::Vector2d heading(std::cos(heading_angle), std::sin(heading_angle));
  const Eigen::Vector2d lateral(-heading.y(), heading.x());
  const Eigen::Vector2d world_velocity = state.segment<2>(3);
  const double linear_velocity = world_velocity.dot(heading);
  const double lateral_velocity = world_velocity.dot(lateral);
  const double half_track_width = config_.track_width / 2.0;
  const double left_slip = control.x() - control.y() * half_track_width -
                           (linear_velocity - state(5) * half_track_width);
  const double right_slip = control.x() + control.y() * half_track_width -
                            (linear_velocity + state(5) * half_track_width);
  const double side_limit = config_.max_traction_acceleration / 2.0;
  const double left_longitudinal =
      std::clamp(config_.linear_response_rate * left_slip / 2.0, -side_limit, side_limit);
  const double right_longitudinal =
      std::clamp(config_.linear_response_rate * right_slip / 2.0, -side_limit, side_limit);
  const double requested_lateral = -config_.lateral_slip_response_rate * lateral_velocity / 2.0;
  const double left_lateral_limit =
      std::sqrt(std::max(0.0, side_limit * side_limit - left_longitudinal * left_longitudinal));
  const double right_lateral_limit =
      std::sqrt(std::max(0.0, side_limit * side_limit - right_longitudinal * right_longitudinal));
  const double lateral_acceleration =
      std::clamp(requested_lateral, -left_lateral_limit, left_lateral_limit) +
      std::clamp(requested_lateral, -right_lateral_limit, right_lateral_limit);
  const double drive_acceleration =
      left_longitudinal + right_longitudinal - config_.linear_drag * linear_velocity;
  const Eigen::Vector2d candidate_world_velocity =
      world_velocity + (heading * drive_acceleration + lateral * lateral_acceleration) * config_.dt;
  const double angular_velocity = state(5);
  const double unconstrained_angular_velocity =
      angular_velocity +
      (config_.angular_response_rate * (right_slip - left_slip) / config_.track_width -
       config_.angular_drag * angular_velocity) *
          config_.dt;
  double next_angular_velocity = std::clamp(
      unconstrained_angular_velocity, -config_.max_angular_velocity, config_.max_angular_velocity);
  double next_heading_angle = WrapAngle(heading_angle + next_angular_velocity * config_.dt);
  Eigen::Vector2d next_heading(std::cos(next_heading_angle), std::sin(next_heading_angle));
  double next_linear_velocity = candidate_world_velocity.dot(next_heading);

  // Wheel-speed bounds also apply to the predicted plant state. Without this
  // projection, momentum can create a twist that is unreachable by the real
  // wheels even if the command itself is feasible.
  for (int iteration = 0; iteration < 2; ++iteration) {
    double left =
        (next_linear_velocity - next_angular_velocity * half_track_width) / config_.wheel_radius;
    double right =
        (next_linear_velocity + next_angular_velocity * half_track_width) / config_.wheel_radius;
    left = std::clamp(left, -config_.max_wheel_velocity, config_.max_wheel_velocity);
    right = std::clamp(right, -config_.max_wheel_velocity, config_.max_wheel_velocity);
    next_linear_velocity = config_.wheel_radius * (left + right) / 2.0;
    next_angular_velocity = config_.wheel_radius * (right - left) / config_.track_width;
    next_heading_angle = WrapAngle(heading_angle + next_angular_velocity * config_.dt);
    next_heading = Eigen::Vector2d(std::cos(next_heading_angle), std::sin(next_heading_angle));
  }
  const Eigen::Vector2d lateral_inertia =
      candidate_world_velocity - next_heading * next_linear_velocity;
  const Eigen::Vector2d next_world_velocity =
      next_heading * next_linear_velocity + config_.lateral_drift_coefficient * lateral_inertia;

  DynamicState next;
  next.segment<2>(0) = state.segment<2>(0) + next_world_velocity * config_.dt;
  next(2) = next_heading_angle;
  next.segment<2>(3) = next_world_velocity;
  next(5) = next_angular_velocity;
  return next;
}

LtvMpc::DynamicState LtvMpc::StateDifference(const DynamicState& left,
                                             const DynamicState& right) const {
  DynamicState difference = left - right;
  difference(2) = WrapAngle(difference(2));
  return difference;
}

void LtvMpc::LinearizeDynamics(const DynamicState& state, const Eigen::Vector2d& control,
                               Eigen::Matrix<double, 6, 6>* state_jacobian,
                               Eigen::Matrix<double, 6, 2>* control_jacobian) const {
  state_jacobian->setZero();
  control_jacobian->setZero();
  const double epsilon = config_.linearization_epsilon;
  for (int index = 0; index < 6; ++index) {
    DynamicState positive = state;
    DynamicState negative = state;
    positive(index) += epsilon;
    negative(index) -= epsilon;
    state_jacobian->col(index) =
        StateDifference(StepModel(positive, control), StepModel(negative, control)) /
        (2.0 * epsilon);
  }
  for (int index = 0; index < 2; ++index) {
    Eigen::Vector2d positive = control;
    Eigen::Vector2d negative = control;
    positive(index) += epsilon;
    negative(index) -= epsilon;
    control_jacobian->col(index) =
        StateDifference(StepModel(state, positive), StepModel(state, negative)) / (2.0 * epsilon);
  }
}

void LtvMpc::BuildQuadratic(const RobotState& state, const BezierTrajectory& reference,
                            double reference_time, Eigen::MatrixXd* hessian,
                            Eigen::VectorXd* gradient) const {
  const int input_count = 2 * config_.horizon_steps;
  *hessian = Eigen::MatrixXd::Zero(input_count, input_count);
  *gradient = Eigen::VectorXd::Zero(input_count);
  const Eigen::Matrix2d input_weight =
      (Eigen::Vector2d(config_.linear_command_weight, config_.angular_command_weight)).asDiagonal();
  const Eigen::Matrix2d delta_input_weight =
      Eigen::Matrix2d::Identity() * config_.control_change_weight;

  const auto reference_state = [&reference](double time) {
    const TrajectoryPoint point = reference.Evaluate(time);
    return DynamicState{point.pose.position.x(), point.pose.position.y(), point.pose.heading,
                        point.velocity.x(),      point.velocity.y(),      point.twist.angular};
  };

  Eigen::Matrix<double, 6, Eigen::Dynamic> sensitivity =
      Eigen::Matrix<double, 6, Eigen::Dynamic>::Zero(6, input_count);
  DynamicState offset = StateDifference(ToDynamicState(state), reference_state(reference_time));
  Eigen::Matrix<double, 2, Eigen::Dynamic> previous_difference =
      Eigen::Matrix<double, 2, Eigen::Dynamic>::Zero(2, input_count);

  for (int index = 0; index < config_.horizon_steps; ++index) {
    const double sample_time = reference_time + static_cast<double>(index) * config_.dt;
    const TrajectoryPoint reference_point = reference.Evaluate(sample_time);
    const DynamicState nominal_state = reference_state(sample_time);
    const Eigen::Vector2d reference_control = ToVector(reference_point.twist);
    const DynamicState next_reference_state = reference_state(sample_time + config_.dt);
    Eigen::Matrix<double, 6, 6> state_jacobian;
    Eigen::Matrix<double, 6, 2> control_jacobian;
    LinearizeDynamics(nominal_state, reference_control, &state_jacobian, &control_jacobian);

    sensitivity = state_jacobian * sensitivity;
    sensitivity.block<6, 2>(0, 2 * index) += control_jacobian;
    offset = state_jacobian * offset - control_jacobian * reference_control +
             StateDifference(StepModel(nominal_state, reference_control), next_reference_state);

    const bool is_terminal = index + 1 == config_.horizon_steps;
    const double heading = next_reference_state(2);
    Eigen::Matrix2d rotation;
    rotation << std::cos(heading), -std::sin(heading), std::sin(heading), std::cos(heading);
    const Eigen::Matrix2d position_weight =
        rotation *
        Eigen::Vector2d(
            is_terminal ? config_.terminal_along_track_weight : config_.along_track_weight,
            is_terminal ? config_.terminal_cross_track_weight : config_.cross_track_weight)
            .asDiagonal() *
        rotation.transpose();
    Eigen::Matrix<double, 6, 6> state_weight = Eigen::Matrix<double, 6, 6>::Zero();
    state_weight.block<2, 2>(0, 0) = position_weight;
    state_weight(2, 2) = is_terminal ? config_.terminal_heading_weight : config_.heading_weight;
    state_weight.block<2, 2>(3, 3) =
        Eigen::Matrix2d::Identity() *
        (is_terminal ? config_.terminal_world_velocity_weight : config_.world_velocity_weight);
    state_weight(5, 5) = is_terminal ? config_.terminal_yaw_rate_weight : config_.yaw_rate_weight;
    *hessian += sensitivity.transpose() * state_weight * sensitivity;
    *gradient += sensitivity.transpose() * state_weight * offset;

    hessian->block<2, 2>(2 * index, 2 * index) += input_weight;
    gradient->segment<2>(2 * index) -= input_weight * reference_control;

    Eigen::Matrix<double, 2, Eigen::Dynamic> difference =
        Eigen::Matrix<double, 2, Eigen::Dynamic>::Zero(2, input_count);
    difference.block<2, 2>(0, 2 * index).setIdentity();
    Eigen::Vector2d difference_offset = Eigen::Vector2d::Zero();
    if (index > 0) {
      difference -= previous_difference;
    } else {
      difference_offset = -ToVector(last_command_);
    }
    *hessian += difference.transpose() * delta_input_weight * difference;
    *gradient += difference.transpose() * delta_input_weight * difference_offset;
    previous_difference = std::move(difference);
  }
}

double LtvMpc::QuadraticCost(const ControlSequence& controls, const Eigen::MatrixXd& hessian,
                             const Eigen::VectorXd& gradient) const {
  Eigen::VectorXd decision(2 * config_.horizon_steps);
  for (int index = 0; index < config_.horizon_steps; ++index) {
    decision.segment<2>(2 * index) = ToVector(controls[index]);
  }
  return decision.dot(hessian * decision) + 2.0 * gradient.dot(decision);
}

}  // namespace differential_drive_mpc
