#ifndef DIFFERENTIAL_DRIVE_MPC_RERUN_VISUALIZER_H_
#define DIFFERENTIAL_DRIVE_MPC_RERUN_VISUALIZER_H_

#include <Eigen/Core>
#include <rerun.hpp>
#include <vector>

#include "differential_drive_mpc/types.h"

namespace differential_drive_mpc {

class RerunVisualizer {
 public:
  explicit RerunVisualizer(double linear_velocity_color_limit);

  void SpawnViewer();
  void LogReference(const std::vector<Eigen::Vector2d>& reference);
  void LogControlPoints(const std::vector<Eigen::Vector2d>& control_points);
  void LogStep(double time, const TrajectoryPoint& target, const RobotState& state, Twist2d command,
               WheelVelocities wheel_velocities, const std::vector<RobotState>& prediction);

 private:
  void ConfigurePlotStyles();

  rerun::RecordingStream recording_;
  std::vector<rerun::Position2D> control_points_;
  std::vector<rerun::Position2D> followed_path_;
  std::vector<double> followed_linear_velocities_;
  std::vector<rerun::Position2D> target_trajectory_;
  double linear_velocity_color_limit_;
};

}  // namespace differential_drive_mpc

#endif  // DIFFERENTIAL_DRIVE_MPC_RERUN_VISUALIZER_H_
