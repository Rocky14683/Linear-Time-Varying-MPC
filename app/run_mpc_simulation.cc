#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "differential_drive_mpc/bezier_trajectory.h"
#include "differential_drive_mpc/differential_drive_robot.h"
#include "differential_drive_mpc/nonlinear_mpc.h"
#include "differential_drive_mpc/rerun_visualizer.h"

int main() {
  using differential_drive_mpc::BezierTrajectory;
  using differential_drive_mpc::DifferentialDriveConfig;
  using differential_drive_mpc::DifferentialDriveRobot;
  using differential_drive_mpc::DisturbanceField;
  using differential_drive_mpc::LtvMpc;
  using differential_drive_mpc::MpcConfig;
  using differential_drive_mpc::Pose2d;
  using differential_drive_mpc::RerunVisualizer;
  using differential_drive_mpc::RobotState;
  using differential_drive_mpc::TrajectoryLimits;
  using differential_drive_mpc::Twist2d;

  const TrajectoryLimits trajectory_limits{.max_linear_velocity = 5.0,
                                           .max_linear_acceleration = 20.0,
                                           .max_angular_velocity = 5.0,
                                           .max_angular_acceleration = 4.0,
                                           .wheel_radius = 0.05,
                                           .track_width = 0.30,
                                           .max_wheel_velocity = 120.0,
                                           .max_wheel_acceleration = 100.0};
  const BezierTrajectory trajectory = BezierTrajectory::Create(
      {Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.5, 0.0), Eigen::Vector2d(2.0, 2.8),
       Eigen::Vector2d(4.8, 2.0), Eigen::Vector2d(0, 2.0), Eigen::Vector2d(-2, 0),
       Eigen::Vector2d(4.8, 1.2), Eigen::Vector2d(8, 1.2)},
      trajectory_limits);

  const double dt = 0.05;
  LtvMpc controller({
      .horizon_steps = 20,
      .dt = dt,
      .optimization_iterations = 50,
      .gradient_step_size = 0.06,
      .max_linear_velocity = trajectory_limits.max_linear_velocity,
      .max_angular_velocity = trajectory_limits.max_angular_velocity,
      .max_linear_acceleration = trajectory_limits.max_linear_acceleration,
      .max_angular_acceleration = trajectory_limits.max_angular_acceleration,
      .wheel_radius = trajectory_limits.wheel_radius,
      .track_width = trajectory_limits.track_width,
      .max_wheel_velocity = trajectory_limits.max_wheel_velocity,
      .max_wheel_acceleration = trajectory_limits.max_wheel_acceleration,
      .linear_response_rate = 20.0,
      .angular_response_rate = 20.0,
      .linear_drag = 0.02,
      .angular_drag = 0.01,
      .lateral_slip_response_rate = 8.0,
      .max_traction_acceleration = 8.0,
      .lateral_drift_coefficient = 0.8,
      .along_track_weight = 10,
      .cross_track_weight = 70,
      .heading_weight = 2,
      .world_velocity_weight = 6,
      .yaw_rate_weight = 0.08,
      .linear_command_weight = 1.8,
      .angular_command_weight = 0.8,
      .control_change_weight = 0.8,
      .terminal_along_track_weight = 100.0,
      .terminal_cross_track_weight = 100.0,
      .terminal_world_velocity_weight = 3.6,
      .terminal_yaw_rate_weight = 0.16,
  });
  DifferentialDriveRobot robot(
      DifferentialDriveConfig{.wheel_radius = 0.05,
                              .track_width = 0.30,
                              .linear_response_rate = 20.0,
                              .angular_response_rate = 20.0,
                              .linear_drag = 0.02,
                              .angular_drag = 0.01,
                              .lateral_slip_response_rate = 8.0,
                              .max_traction_acceleration = 8.0,
                              .lateral_drift_coefficient = 0.8},
      RobotState{.pose = Pose2d{.position = Eigen::Vector2d(0.0, -0.25), .heading = 0.0}},
      DisturbanceField{.base_linear_friction = 0.03,
                       .x_friction_gradient = 0.045,
                       .y_friction_gradient = 0.010,
                       .external_acceleration = Eigen::Vector2d(0.0, 0.0),
                       .external_angular_acceleration = 0.2});

  RerunVisualizer visualizer(trajectory_limits.max_linear_velocity);
  visualizer.SpawnViewer();
  visualizer.LogReference(trajectory.SamplePositionsByDistance(300));
  visualizer.LogControlPoints(trajectory.control_points());

  const int steps = static_cast<int>(std::ceil(trajectory.duration() / dt));
  for (int step = 0; step <= steps; ++step) {
    const double time = std::min(step * dt, trajectory.duration());
    const Twist2d command = controller.Solve(robot.state(), trajectory, time);
    robot.SetCommand(command);
    robot.Step(dt);
    visualizer.LogStep(time, trajectory.Evaluate(time), robot.state(), command,
                       robot.commanded_wheel_velocities(), controller.predicted_states());
  }

  const Eigen::Vector2d goal = trajectory.Evaluate(trajectory.duration()).pose.position;
  const std::vector<differential_drive_mpc::TrajectoryPoint> reference_samples =
      trajectory.Sample(10001);
  double max_reference_linear_velocity = 0.0;
  double max_reference_linear_acceleration = 0.0;
  double max_reference_angular_velocity = 0.0;
  double max_reference_angular_acceleration = 0.0;
  for (const differential_drive_mpc::TrajectoryPoint& point : reference_samples) {
    max_reference_linear_velocity =
        std::max(max_reference_linear_velocity, std::abs(point.twist.linear));
    max_reference_linear_acceleration =
        std::max(max_reference_linear_acceleration, std::abs(point.acceleration_twist.linear));
    max_reference_angular_velocity =
        std::max(max_reference_angular_velocity, std::abs(point.twist.angular));
    max_reference_angular_acceleration =
        std::max(max_reference_angular_acceleration, std::abs(point.acceleration_twist.angular));
  }
  std::cout << "Simulation complete. Final position error: "
            << (robot.state().pose.position - goal).norm() << " m\n";
  std::cout << "Final velocity: " << robot.state().linear_velocity << " m/s, "
            << robot.state().angular_velocity << " rad/s\n";
  std::cout << "Reference extrema (v, a, omega, alpha): " << max_reference_linear_velocity << ", "
            << max_reference_linear_acceleration << ", " << max_reference_angular_velocity << ", "
            << max_reference_angular_acceleration << "\n";
  return 0;
}
