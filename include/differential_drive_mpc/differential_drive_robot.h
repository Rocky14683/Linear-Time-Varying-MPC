#ifndef DIFFERENTIAL_DRIVE_MPC_DIFFERENTIAL_DRIVE_ROBOT_H_
#define DIFFERENTIAL_DRIVE_MPC_DIFFERENTIAL_DRIVE_ROBOT_H_

#include <Eigen/Core>

#include "differential_drive_mpc/types.h"

namespace differential_drive_mpc {

// Parameters for the command-to-motion plant model.
struct DifferentialDriveConfig {
  double wheel_radius = 0.05;
  double track_width = 0.30;
  double linear_response_rate = 12.0;
  double angular_response_rate = 14.0;
  double linear_drag = 0.0;
  double angular_drag = 0.0;
  double lateral_slip_response_rate = 0.0;
  double max_traction_acceleration = 50.0;
  // Fraction of lateral world-velocity inertia retained per simulation step.
  // 0.0 enforces ideal no-slip motion; 1.0 retains all lateral inertia.
  double lateral_drift_coefficient = 0.0;
};

// A deterministic disturbance field. Friction changes with world position;
// external_acceleration is a world-frame acceleration applied to the chassis.
struct DisturbanceField {
  double base_linear_friction = 0.0;
  double x_friction_gradient = 0.0;
  double y_friction_gradient = 0.0;
  double max_linear_friction = 2.0;
  Eigen::Vector2d external_acceleration = Eigen::Vector2d::Zero();
  double external_angular_acceleration = 0.0;
};

class DifferentialDriveRobot {
 public:
  DifferentialDriveRobot(DifferentialDriveConfig config, RobotState initial_state,
                         DisturbanceField disturbances = {});

  void SetCommand(Twist2d command);
  void Step(double dt);

  WheelVelocities CommandToWheelVelocities(Twist2d command) const;
  Twist2d WheelVelocitiesToCommand(WheelVelocities wheel_velocities) const;

  const RobotState& state() const { return state_; }
  const Twist2d& command() const { return command_; }
  WheelVelocities commanded_wheel_velocities() const;
  double last_linear_friction() const { return last_linear_friction_; }

 private:
  double LinearFrictionAt(const Eigen::Vector2d& position) const;

  DifferentialDriveConfig config_;
  DisturbanceField disturbances_;
  RobotState state_;
  Twist2d command_;
  double last_linear_friction_ = 0.0;
};

}  // namespace differential_drive_mpc

#endif  // DIFFERENTIAL_DRIVE_MPC_DIFFERENTIAL_DRIVE_ROBOT_H_
