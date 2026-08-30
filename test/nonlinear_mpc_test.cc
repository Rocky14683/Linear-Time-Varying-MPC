#include "differential_drive_mpc/nonlinear_mpc.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <cmath>

#include "differential_drive_mpc/bezier_trajectory.h"
#include "differential_drive_mpc/differential_drive_robot.h"

namespace differential_drive_mpc {
namespace {

BezierTrajectory MakeStraightTrajectory() {
  return BezierTrajectory::Create({Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.0, 0.0),
                                   Eigen::Vector2d(3.0, 0.0), Eigen::Vector2d(4.0, 0.0)},
                                  {});
}

TEST(AngleErrorTest, UsesTheShortestSignedRotationAcrossThePiBoundary) {
  constexpr double kSmallAngle = 0.02;
  EXPECT_NEAR(ShortestAngularDifference(-kPi + kSmallAngle, kPi - kSmallAngle),
              2.0 * kSmallAngle, 1e-12);
  EXPECT_NEAR(ShortestAngularDifference(kPi - kSmallAngle, -kPi + kSmallAngle),
              -2.0 * kSmallAngle, 1e-12);
}

TEST(LtvMpcTest, ProducesAFeasibleForwardCommandAndPrediction) {
  MpcConfig config{.horizon_steps = 10,
                   .dt = 0.1,
                   .optimization_iterations = 10,
                   .gradient_step_size = 0.08,
                   .max_linear_acceleration = 2.0,
                   .max_angular_acceleration = 3.0};
  LtvMpc controller(config);
  const BezierTrajectory trajectory = MakeStraightTrajectory();

  const Twist2d command = controller.Solve({}, trajectory, 0.0);
  EXPECT_GT(command.linear, 0.0);
  EXPECT_LE(command.linear, config.max_linear_acceleration * config.dt + 1e-8);
  EXPECT_LE(std::abs(command.angular), config.max_angular_acceleration * config.dt + 1e-8);
  ASSERT_EQ(controller.predicted_states().size(), config.horizon_steps);
  EXPECT_GT(controller.predicted_states().back().pose.position.x(), 0.0);
  EXPECT_TRUE(std::isfinite(controller.last_cost()));
}

TEST(LtvMpcTest, EnforcesWheelRateBoundsAcrossThePredictionHorizon) {
  MpcConfig config{.horizon_steps = 12,
                   .dt = 0.05,
                   .optimization_iterations = 30,
                   .gradient_step_size = 0.06,
                   .max_linear_velocity = 10.0,
                   .max_angular_velocity = 10.0,
                   .max_linear_acceleration = 100.0,
                   .max_angular_acceleration = 100.0,
                   .wheel_radius = 0.10,
                   .track_width = 0.40,
                   .max_wheel_velocity = 3.0,
                   .max_wheel_acceleration = 100.0};
  LtvMpc controller(config);
  const BezierTrajectory trajectory =
      BezierTrajectory::Create({Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(0.8, 0.0),
                                Eigen::Vector2d(0.0, 0.8), Eigen::Vector2d(0.8, 0.8)},
                               {.wheel_radius = config.wheel_radius,
                                .track_width = config.track_width,
                                .max_wheel_velocity = config.max_wheel_velocity,
                                .max_wheel_acceleration = config.max_wheel_acceleration});

  const Twist2d command = controller.Solve({}, trajectory, trajectory.duration() * 0.5);
  const auto check_wheels = [&config](Twist2d twist) {
    const double half_track = config.track_width / 2.0;
    EXPECT_LE(std::abs((twist.linear - twist.angular * half_track) / config.wheel_radius),
              config.max_wheel_velocity + 1e-9);
    EXPECT_LE(std::abs((twist.linear + twist.angular * half_track) / config.wheel_radius),
              config.max_wheel_velocity + 1e-9);
  };
  check_wheels(command);
  for (const RobotState& prediction : controller.predicted_states()) {
    check_wheels({.linear = prediction.linear_velocity, .angular = prediction.angular_velocity});
  }
}

TEST(LtvMpcTest, CompensatesCrossTrackErrorWithLateralInertia) {
  MpcConfig config{.horizon_steps = 12,
                   .dt = 0.05,
                   .optimization_iterations = 12,
                   .gradient_step_size = 0.06,
                   .max_linear_acceleration = 3.0,
                   .max_angular_acceleration = 5.0,
                   .linear_response_rate = 20.0,
                   .angular_response_rate = 20.0,
                   .lateral_drift_coefficient = 0.7,
                   .along_track_weight = 5.0,
                   .cross_track_weight = 50.0,
                   .terminal_along_track_weight = 20.0,
                   .terminal_cross_track_weight = 100.0};
  LtvMpc controller(config);
  const BezierTrajectory trajectory = MakeStraightTrajectory();
  DifferentialDriveRobot robot({.linear_response_rate = 20.0,
                                .angular_response_rate = 20.0,
                                .lateral_drift_coefficient = 0.7},
                               {.pose = {.position = Eigen::Vector2d(0.0, 0.5), .heading = 0.0}});

  const double initial_error = std::abs(robot.state().pose.position.y());
  for (int step = 0; step < 30; ++step) {
    const Twist2d command = controller.Solve(robot.state(), trajectory, step * config.dt);
    robot.SetCommand(command);
    robot.Step(config.dt);
  }
  EXPECT_LT(std::abs(robot.state().pose.position.y()), initial_error);
}

}  // namespace
}  // namespace differential_drive_mpc
