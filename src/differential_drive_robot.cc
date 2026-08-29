#include "differential_drive_mpc/differential_drive_robot.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace differential_drive_mpc {

DifferentialDriveRobot::DifferentialDriveRobot(DifferentialDriveConfig config,
                                               RobotState initial_state,
                                               DisturbanceField disturbances)
    : config_(config), disturbances_(disturbances), state_(initial_state) {
  if (config_.wheel_radius <= 0.0 || config_.track_width <= 0.0 ||
      config_.linear_response_rate < 0.0 || config_.angular_response_rate < 0.0) {
    throw std::invalid_argument("Wheel geometry must be positive and response rates nonnegative.");
  }
  if (!std::isfinite(config_.lateral_drift_coefficient) ||
      config_.lateral_drift_coefficient < 0.0 || config_.lateral_drift_coefficient > 1.0) {
    throw std::invalid_argument("Lateral drift coefficient must be in [0, 1].");
  }
  if (config_.lateral_slip_response_rate < 0.0 || config_.max_traction_acceleration <= 0.0) {
    throw std::invalid_argument("Slip response must be nonnegative and traction positive.");
  }
  const Eigen::Vector2d initial_heading(std::cos(state_.pose.heading),
                                        std::sin(state_.pose.heading));
  if (state_.world_velocity.isZero() && state_.linear_velocity != 0.0) {
    state_.world_velocity = initial_heading * state_.linear_velocity;
  } else {
    state_.linear_velocity = state_.world_velocity.dot(initial_heading);
  }
}

void DifferentialDriveRobot::SetCommand(Twist2d command) { command_ = command; }

void DifferentialDriveRobot::Step(double dt) {
  if (dt <= 0.0) {
    throw std::invalid_argument("Simulation time step must be positive.");
  }

  last_linear_friction_ = LinearFrictionAt(state_.pose.position);
  const Eigen::Vector2d heading(std::cos(state_.pose.heading), std::sin(state_.pose.heading));
  const Eigen::Vector2d lateral(-heading.y(), heading.x());
  state_.linear_velocity = state_.world_velocity.dot(heading);
  const double lateral_velocity = state_.world_velocity.dot(lateral);
  const double half_track = config_.track_width / 2.0;
  const double left_slip = command_.linear - command_.angular * half_track -
                           (state_.linear_velocity - state_.angular_velocity * half_track);
  const double right_slip = command_.linear + command_.angular * half_track -
                            (state_.linear_velocity + state_.angular_velocity * half_track);
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
      left_longitudinal + right_longitudinal -
      (config_.linear_drag + last_linear_friction_) * state_.linear_velocity;
  const double angular_acceleration =
      config_.angular_response_rate * (right_slip - left_slip) / config_.track_width -
      config_.angular_drag * state_.angular_velocity + disturbances_.external_angular_acceleration;

  // World velocity is the translational state. Wheel drive contributes along
  // the current heading, while external acceleration acts directly in world
  // coordinates. This candidate retains momentum from prior time steps.
  const Eigen::Vector2d candidate_world_velocity =
      state_.world_velocity + (heading * drive_acceleration + lateral * lateral_acceleration +
                               disturbances_.external_acceleration) *
                                  dt;
  state_.angular_velocity += angular_acceleration * dt;
  state_.pose.heading = WrapAngle(state_.pose.heading + state_.angular_velocity * dt);
  const Eigen::Vector2d updated_heading(std::cos(state_.pose.heading),
                                        std::sin(state_.pose.heading));
  const double projected_linear_velocity = candidate_world_velocity.dot(updated_heading);
  const Eigen::Vector2d lateral_inertia =
      candidate_world_velocity - updated_heading * projected_linear_velocity;
  state_.linear_velocity = projected_linear_velocity;
  state_.world_velocity = updated_heading * projected_linear_velocity +
                          config_.lateral_drift_coefficient * lateral_inertia;
  state_.pose.position += state_.world_velocity * dt;
}

WheelVelocities DifferentialDriveRobot::CommandToWheelVelocities(Twist2d command) const {
  return {
      .left = (command.linear - command.angular * config_.track_width / 2.0) / config_.wheel_radius,
      .right =
          (command.linear + command.angular * config_.track_width / 2.0) / config_.wheel_radius};
}

Twist2d DifferentialDriveRobot::WheelVelocitiesToCommand(WheelVelocities wheel_velocities) const {
  return {.linear = config_.wheel_radius * (wheel_velocities.left + wheel_velocities.right) / 2.0,
          .angular = config_.wheel_radius * (wheel_velocities.right - wheel_velocities.left) /
                     config_.track_width};
}

WheelVelocities DifferentialDriveRobot::commanded_wheel_velocities() const {
  return CommandToWheelVelocities(command_);
}

double DifferentialDriveRobot::LinearFrictionAt(const Eigen::Vector2d& position) const {
  const double friction = disturbances_.base_linear_friction +
                          disturbances_.x_friction_gradient * position.x() +
                          disturbances_.y_friction_gradient * position.y();
  return std::clamp(friction, 0.0, disturbances_.max_linear_friction);
}

}  // namespace differential_drive_mpc
