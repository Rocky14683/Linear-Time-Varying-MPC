#ifndef DIFFERENTIAL_DRIVE_MPC_TYPES_H_
#define DIFFERENTIAL_DRIVE_MPC_TYPES_H_

#include <Eigen/Core>
#include <cmath>

namespace differential_drive_mpc {

constexpr double kPi = 3.14159265358979323846;

inline double WrapAngle(double angle_radians) { return std::remainder(angle_radians, 2.0 * kPi); }

struct Pose2d {
  Eigen::Vector2d position = Eigen::Vector2d::Zero();
  double heading = 0.0;
};

struct Twist2d {
  double linear = 0.0;
  double angular = 0.0;
};

struct WheelVelocities {
  double left = 0.0;
  double right = 0.0;
};

struct RobotState {
  Pose2d pose;
  // Chassis translational velocity resolved in world coordinates. It is the
  // velocity integrated by the plant and may differ from the no-slip heading
  // direction when the drift model is enabled.
  Eigen::Vector2d world_velocity = Eigen::Vector2d::Zero();
  double linear_velocity = 0.0;
  double angular_velocity = 0.0;
};

struct TrajectoryPoint {
  double time = 0.0;
  Pose2d pose;
  Eigen::Vector2d velocity = Eigen::Vector2d::Zero();
  Eigen::Vector2d acceleration = Eigen::Vector2d::Zero();
  Twist2d twist;
  Twist2d acceleration_twist;
};

}  // namespace differential_drive_mpc

#endif  // DIFFERENTIAL_DRIVE_MPC_TYPES_H_
