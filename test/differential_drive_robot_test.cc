#include "differential_drive_mpc/differential_drive_robot.h"

#include <gtest/gtest.h>

namespace differential_drive_mpc {
namespace {

TEST(DifferentialDriveRobotTest, ConvertsCommandsToWheelVelocitiesAndBack) {
  DifferentialDriveRobot robot({.wheel_radius = 0.1, .track_width = 0.4}, {});
  const Twist2d command{.linear = 1.2, .angular = -0.5};

  const WheelVelocities wheels = robot.CommandToWheelVelocities(command);
  EXPECT_NEAR(wheels.left, 13.0, 1e-12);
  EXPECT_NEAR(wheels.right, 11.0, 1e-12);
  const Twist2d reconstructed = robot.WheelVelocitiesToCommand(wheels);
  EXPECT_NEAR(reconstructed.linear, command.linear, 1e-12);
  EXPECT_NEAR(reconstructed.angular, command.angular, 1e-12);
}

TEST(DifferentialDriveRobotTest, IntegratesPoseFromVelocityCommand) {
  DifferentialDriveRobot robot({.linear_response_rate = 10.0, .angular_response_rate = 10.0}, {});
  robot.SetCommand({.linear = 1.0, .angular = 0.0});
  robot.Step(0.1);

  EXPECT_NEAR(robot.state().linear_velocity, 1.0, 1e-12);
  EXPECT_NEAR(robot.state().pose.position.x(), 0.1, 1e-12);
  EXPECT_NEAR(robot.state().pose.position.y(), 0.0, 1e-12);
}

TEST(DifferentialDriveRobotTest, ProducesAndIntegratesWorldFrameDriftVelocity) {
  const DifferentialDriveConfig config{.linear_response_rate = 10.0,
                                       .angular_response_rate = 10.0,
                                       .lateral_drift_coefficient = 0.2};
  DifferentialDriveRobot robot(config, {.linear_velocity = 1.0});
  EXPECT_NEAR(robot.state().world_velocity.x(), 1.0, 1e-12);
  EXPECT_NEAR(robot.state().world_velocity.y(), 0.0, 1e-12);
  robot.SetCommand({.linear = 1.0, .angular = 1.0});
  robot.Step(0.1);

  const double heading = 0.1;
  const Eigen::Vector2d forward(std::cos(heading), std::sin(heading));
  const Eigen::Vector2d candidate_world_velocity(1.0, 0.0);
  const Eigen::Vector2d lateral_inertia =
      candidate_world_velocity - forward * candidate_world_velocity.dot(forward);
  const Eigen::Vector2d expected_world_velocity =
      forward * candidate_world_velocity.dot(forward) + 0.2 * lateral_inertia;
  const Eigen::Vector2d expected_position = expected_world_velocity * 0.1;
  EXPECT_NEAR(robot.state().world_velocity.x(), expected_world_velocity.x(), 1e-12);
  EXPECT_NEAR(robot.state().world_velocity.y(), expected_world_velocity.y(), 1e-12);
  EXPECT_NEAR(robot.state().pose.position.x(), expected_position.x(), 1e-12);
  EXPECT_NEAR(robot.state().pose.position.y(), expected_position.y(), 1e-12);

  // The lateral component persists into the next step instead of being an
  // instantaneous, algebraic offset.
  robot.Step(0.1);
  EXPECT_NE(robot.state().world_velocity.dot(Eigen::Vector2d(-std::sin(robot.state().pose.heading),
                                                             std::cos(robot.state().pose.heading))),
            0.0);
}

TEST(DifferentialDriveRobotTest, AppliesPositionDependentFriction) {
  const DifferentialDriveConfig config{.linear_response_rate = 0.0, .angular_response_rate = 0.0};
  const DisturbanceField disturbances{.base_linear_friction = 0.2, .x_friction_gradient = 0.5};
  const RobotState at_origin{.linear_velocity = 1.0};
  const RobotState at_positive_x{.pose = {.position = Eigen::Vector2d(1.0, 0.0)},
                                 .linear_velocity = 1.0};
  DifferentialDriveRobot low_friction(config, at_origin, disturbances);
  DifferentialDriveRobot high_friction(config, at_positive_x, disturbances);

  low_friction.Step(0.1);
  high_friction.Step(0.1);

  EXPECT_LT(high_friction.state().linear_velocity, low_friction.state().linear_velocity);
  EXPECT_GT(high_friction.last_linear_friction(), low_friction.last_linear_friction());
}

}  // namespace
}  // namespace differential_drive_mpc
