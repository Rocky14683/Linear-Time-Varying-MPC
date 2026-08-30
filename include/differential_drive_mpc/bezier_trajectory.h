#ifndef DIFFERENTIAL_DRIVE_MPC_BEZIER_TRAJECTORY_H_
#define DIFFERENTIAL_DRIVE_MPC_BEZIER_TRAJECTORY_H_

#include <Eigen/Core>
#include <vector>

#include "differential_drive_mpc/types.h"

namespace differential_drive_mpc {

// Type-erased reference consumed by the MPC. Concrete trajectory types can
// provide their own path geometry while sharing the controller.
class Trajectory {
 public:
  virtual ~Trajectory() = default;
  virtual TrajectoryPoint Evaluate(double time) const = 0;
  virtual double duration() const = 0;
};

// Static adapter for concrete trajectory implementations. The virtual surface
// is intentionally limited to the controller-facing reference; path-specific
// work remains statically dispatched through Derived.
template <typename Derived>
class TrajectoryBase : public Trajectory {
 public:
  TrajectoryPoint Evaluate(double time) const final {
    return static_cast<const Derived*>(this)->EvaluateImpl(time);
  }
  double duration() const final { return static_cast<const Derived*>(this)->duration_impl(); }
};

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

// A C2 interpolating cubic B-spline through the supplied waypoints, with a
// forward/backward, minimum-time scaling pass that respects chassis and wheel
// limits.
class CubicBSplineTrajectory : public TrajectoryBase<CubicBSplineTrajectory> {
 public:
  // Creates the minimum-time trajectory through at least four waypoints.
  static CubicBSplineTrajectory Create(std::vector<Eigen::Vector2d> waypoints,
                                       const TrajectoryLimits& limits, bool reverse = false);

  TrajectoryPoint EvaluateImpl(double time) const;
  std::vector<TrajectoryPoint> Sample(int sample_count) const;
  // Returns spatially uniform samples for rendering the geometric path. Unlike
  // Sample(), these positions are uniformly spaced by arc length, not time.
  std::vector<Eigen::Vector2d> SamplePositionsByDistance(int sample_count) const;

  const std::vector<Eigen::Vector2d>& control_points() const { return waypoints_; }
  double duration_impl() const { return duration_; }
  double total_length() const { return total_length_; }
  const TrajectoryLimits& limits() const { return limits_; }
  bool reverse() const { return reverse_; }

 private:
  struct ArcLengthSample {
    double path_parameter = 0.0;
    double distance = 0.0;
  };

  struct TimeProfileSample {
    double distance = 0.0;
    double time = 0.0;
    double speed = 0.0;
  };

  struct ProgressState {
    double distance = 0.0;
    double speed = 0.0;
    double acceleration = 0.0;
  };

  CubicBSplineTrajectory(std::vector<Eigen::Vector2d> waypoints, TrajectoryLimits limits,
                          bool reverse);

  void BuildInterpolatingSpline();
  void BuildArcLengthTable();
  void BuildMinimumTimeProfile();
  Eigen::Vector2d EvaluateDerivative(double path_parameter, int order) const;
  double ParameterAtDistance(double distance) const;
  double CurvatureAtDistance(double distance) const;
  double CurvatureDerivativeAtDistance(double distance) const;
  double VelocityLimitAtDistance(double distance) const;
  bool AccelerationBoundsAtDistance(double distance, double speed, double* lower,
                                    double* upper) const;
  ProgressState ProgressAtTime(double time) const;

  std::vector<ArcLengthSample> arc_length_table_;
  std::vector<TimeProfileSample> time_profile_;
  std::vector<Eigen::Vector2d> waypoints_;
  std::vector<Eigen::Vector2d> spline_control_points_;
  std::vector<double> knots_;
  double duration_;
  double total_length_ = 0.0;
  TrajectoryLimits limits_;
  bool reverse_ = false;
};

// Transitional source compatibility for existing callers. New code should use
// CubicBSplineTrajectory explicitly.
using BezierTrajectory = CubicBSplineTrajectory;

}  // namespace differential_drive_mpc

#endif  // DIFFERENTIAL_DRIVE_MPC_BEZIER_TRAJECTORY_H_
