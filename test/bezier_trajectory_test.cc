#include "differential_drive_mpc/bezier_trajectory.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <algorithm>
#include <vector>

namespace differential_drive_mpc {
namespace {

TEST(BezierTrajectoryTest, ScalesDurationAndRespectsAllMotionLimits) {
  const TrajectoryLimits limits{.max_linear_velocity = 1.0,
                                .max_linear_acceleration = 1.0,
                                .max_angular_velocity = 1.0,
                                .max_angular_acceleration = 2.0};
  const BezierTrajectory trajectory =
      BezierTrajectory::Create({Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.0, 0.0),
                                Eigen::Vector2d(1.0, 2.0), Eigen::Vector2d(3.0, 2.0)},
                               1.0, limits);

  EXPECT_GT(trajectory.duration(), 1.0);
  for (const TrajectoryPoint& point : trajectory.Sample(400)) {
    EXPECT_LE(point.twist.linear, limits.max_linear_velocity + 1e-9);
    EXPECT_LE(std::abs(point.acceleration_twist.linear), limits.max_linear_acceleration + 1e-9);
    EXPECT_LE(std::abs(point.twist.angular), limits.max_angular_velocity + 1e-9);
    EXPECT_LE(std::abs(point.acceleration_twist.angular), limits.max_angular_acceleration + 1e-9);
  }
}

TEST(BezierTrajectoryTest, CapsComposedWheelVelocityOnTightCurves) {
  // This deliberately leaves the independent v and omega limits loose. The
  // wheel bound must still slow the time profile whenever turning combines
  // translational and rotational motion at one wheel.
  const TrajectoryLimits limits{.max_linear_velocity = 20.0,
                                .max_linear_acceleration = 100.0,
                                .max_angular_velocity = 20.0,
                                .max_angular_acceleration = 100.0,
                                .wheel_radius = 0.10,
                                .track_width = 0.40,
                                .max_wheel_velocity = 3.0,
                                .max_wheel_acceleration = 100.0};
  const BezierTrajectory trajectory =
      BezierTrajectory::Create({Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(0.8, 0.0),
                                Eigen::Vector2d(0.0, 0.8), Eigen::Vector2d(0.8, 0.8)},
                               1.0, limits);

  for (const TrajectoryPoint& point : trajectory.Sample(1000)) {
    const double half_track = limits.track_width / 2.0;
    const double left =
        (point.twist.linear - point.twist.angular * half_track) / limits.wheel_radius;
    const double right =
        (point.twist.linear + point.twist.angular * half_track) / limits.wheel_radius;
    EXPECT_LE(std::abs(left), limits.max_wheel_velocity + 1e-9);
    EXPECT_LE(std::abs(right), limits.max_wheel_velocity + 1e-9);
  }
}

TEST(BezierTrajectoryTest, GeneratesAHeadingTangentToThePath) {
  const BezierTrajectory trajectory =
      BezierTrajectory::Create({Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.0, 0.0),
                                Eigen::Vector2d(2.0, 1.0), Eigen::Vector2d(3.0, 1.0)},
                               4.0, {});

  const TrajectoryPoint point = trajectory.Evaluate(trajectory.duration() * 0.4);
  const Eigen::Vector2d lateral(-std::sin(point.pose.heading), std::cos(point.pose.heading));
  EXPECT_NEAR(point.velocity.dot(lateral), 0.0, 1e-10);
  EXPECT_NEAR(trajectory.Evaluate(0.0).pose.position.x(), 0.0, 1e-12);
  EXPECT_NEAR(trajectory.Evaluate(trajectory.duration()).pose.position.y(), 1.0, 1e-12);
}

TEST(BezierTrajectoryTest, UsesArcLengthForItsLinearSpeedProfile) {
  const BezierTrajectory trajectory =
      BezierTrajectory::Create({Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(0.2, 0.0),
                                Eigen::Vector2d(4.0, 3.0), Eigen::Vector2d(5.0, 3.0)},
                               8.0, {});

  const TrajectoryPoint midpoint = trajectory.Evaluate(trajectory.duration() * 0.5);
  // The minimum-jerk distance profile has q'(0.5) = 1.875. Consequently,
  // the physical speed is independent of the Bezier control parameter rate.
  const double expected_speed = 1.875 * trajectory.total_length() / trajectory.duration();
  EXPECT_NEAR(midpoint.twist.linear, expected_speed, 1e-8);
}

TEST(BezierTrajectoryTest, SamplesDisplayedPathUniformlyByDistance) {
  const BezierTrajectory trajectory =
      BezierTrajectory::Create({Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(0.1, 0.0),
                                Eigen::Vector2d(8.0, 0.0), Eigen::Vector2d(10.0, 0.0)},
                               8.0, {});

  const std::vector<Eigen::Vector2d> positions = trajectory.SamplePositionsByDistance(11);
  ASSERT_EQ(positions.size(), 11U);
  for (size_t index = 0; index < positions.size(); ++index) {
    EXPECT_NEAR(positions[index].x(), static_cast<double>(index), 1e-8);
    EXPECT_NEAR(positions[index].y(), 0.0, 1e-12);
  }
}

TEST(BezierTrajectoryTest, StartsAndEndsAtRestWithoutViolatingAccelerationLimits) {
  const TrajectoryLimits limits{.max_linear_velocity = 1.0,
                                .max_linear_acceleration = 0.8,
                                .max_angular_velocity = 1.0,
                                .max_angular_acceleration = 1.5};
  const BezierTrajectory trajectory =
      BezierTrajectory::Create({Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.0, 0.0),
                                Eigen::Vector2d(1.0, 2.0), Eigen::Vector2d(3.0, 2.0)},
                               1.0, limits);

  const TrajectoryPoint start = trajectory.Evaluate(0.0);
  const TrajectoryPoint end = trajectory.Evaluate(trajectory.duration());
  EXPECT_NEAR(start.twist.linear, 0.0, 1e-12);
  EXPECT_NEAR(start.twist.angular, 0.0, 1e-12);
  EXPECT_NEAR(end.twist.linear, 0.0, 1e-12);
  EXPECT_NEAR(end.twist.angular, 0.0, 1e-12);
  EXPECT_NEAR(end.acceleration_twist.linear, 0.0, 1e-12);
  EXPECT_NEAR(end.acceleration_twist.angular, 0.0, 1e-12);

  bool has_nonzero_speed = false;
  for (const TrajectoryPoint& point : trajectory.Sample(400)) {
    has_nonzero_speed = has_nonzero_speed || point.twist.linear > 0.1;
    EXPECT_LE(std::abs(point.acceleration_twist.linear), limits.max_linear_acceleration + 1e-9);
    EXPECT_LE(std::abs(point.acceleration_twist.angular), limits.max_angular_acceleration + 1e-9);
  }
  EXPECT_TRUE(has_nonzero_speed);
}

TEST(BezierTrajectoryTest, RejectsACusp) {
  EXPECT_THROW(BezierTrajectory::Create({Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(0.0, 0.0),
                                         Eigen::Vector2d(1.0, 0.0), Eigen::Vector2d(2.0, 0.0)},
                                        1.0, {}),
               std::invalid_argument);
}

}  // namespace
}  // namespace differential_drive_mpc
