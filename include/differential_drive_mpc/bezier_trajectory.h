#ifndef DIFFERENTIAL_DRIVE_MPC_BEZIER_TRAJECTORY_H_
#define DIFFERENTIAL_DRIVE_MPC_BEZIER_TRAJECTORY_H_

#include <Eigen/Core>
#include <vector>

#include "differential_drive_mpc/types.h"

namespace differential_drive_mpc {

// Limits for a time-parameterized, differentially feasible reference path.
struct TrajectoryLimits {
  double max_linear_velocity = 10.0;
  double max_linear_acceleration = 5.0;
  double max_angular_velocity = 2.0;
  double max_angular_acceleration = 4.0;
  // Differential-drive geometry and actuator limits.  A valid twist must
  // satisfy |(v +/- omega * track_width / 2) / wheel_radius| <=
  // max_wheel_velocity for both wheels; limiting v and omega independently
  // does not guarantee this.
  double wheel_radius = 0.05;
  double track_width = 0.30;
  double max_wheel_velocity = 200.0;
  double max_wheel_acceleration = 1000.0;
};

// A Bezier path reparameterized by arc length, with a duration selected to
// respect the supplied motion limits.
class BezierTrajectory {
 public:
  // Creates a trajectory from at least four control points. The duration is
  // increased when needed to satisfy all limits.
  static BezierTrajectory Create(std::vector<Eigen::Vector2d> control_points,
                                 double requested_duration, const TrajectoryLimits& limits);

  TrajectoryPoint Evaluate(double time) const;
  std::vector<TrajectoryPoint> Sample(int sample_count) const;
  // Returns spatially uniform samples for rendering the geometric path. Unlike
  // Sample(), these positions are uniformly spaced by arc length, not time.
  std::vector<Eigen::Vector2d> SamplePositionsByDistance(int sample_count) const;

  const std::vector<Eigen::Vector2d>& control_points() const { return control_points_; }
  double duration() const { return duration_; }
  double total_length() const { return total_length_; }
  const TrajectoryLimits& limits() const { return limits_; }

 private:
  struct ArcLengthSample {
    double bezier_parameter = 0.0;
    double distance = 0.0;
  };

  BezierTrajectory(std::vector<Eigen::Vector2d> control_points, double duration,
                   TrajectoryLimits limits);

  void BuildArcLengthTable();
  Eigen::Vector2d EvaluateDerivative(double bezier_parameter, int order) const;
  double ParameterAtDistance(double distance) const;
  double FindRequiredTimeScale() const;

  std::vector<ArcLengthSample> arc_length_table_;
  std::vector<Eigen::Vector2d> control_points_;
  double duration_;
  double total_length_ = 0.0;
  TrajectoryLimits limits_;
};

}  // namespace differential_drive_mpc

#endif  // DIFFERENTIAL_DRIVE_MPC_BEZIER_TRAJECTORY_H_
