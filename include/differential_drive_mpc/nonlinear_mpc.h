#ifndef DIFFERENTIAL_DRIVE_MPC_NONLINEAR_MPC_H_
#define DIFFERENTIAL_DRIVE_MPC_NONLINEAR_MPC_H_

#include <vector>

#include "differential_drive_mpc/bezier_trajectory.h"
#include "differential_drive_mpc/types.h"

namespace differential_drive_mpc {

struct MpcConfig {
  int horizon_steps = 20;
  double dt = 0.05;
  // Projected-gradient iterations for the condensed constrained QP.
  int optimization_iterations = 18;
  double gradient_step_size = 0.10;
  double max_linear_velocity = 2.0;
  double max_angular_velocity = 2.5;
  double max_linear_acceleration = 2.5;
  double max_angular_acceleration = 5.0;
  double wheel_radius = 0.05;
  double track_width = 0.30;
  double max_wheel_velocity = 40.0;
  double max_wheel_acceleration = 200.0;
  // Prediction-model parameters. Match these to DifferentialDriveConfig so
  // the MPC can anticipate the plant's persistent world-frame drift.
  double linear_response_rate = 12.0;
  double angular_response_rate = 14.0;
  double linear_drag = 0.0;
  double angular_drag = 0.0;
  double lateral_slip_response_rate = 0.0;
  double max_traction_acceleration = 50.0;
  double lateral_drift_coefficient = 0.0;
  double linearization_epsilon = 1e-5;
  // Position error weights in the reference-path frame. The cross-track term
  // is deliberately higher so the robot prioritizes staying on the curve.
  double along_track_weight = 10.0;
  double cross_track_weight = 50.0;
  double heading_weight = 3.0;
  // Dynamic-state weights used to compensate velocity drift.
  double world_velocity_weight = 4.0;
  double yaw_rate_weight = 1.0;
  // Input-deviation and slew penalties keep the correction close to the
  // wheel-feasible feedforward command and avoid actuator chatter.
  double linear_command_weight = 1.0;
  double angular_command_weight = 0.5;
  double control_change_weight = 1;
  double terminal_along_track_weight = 20.0;
  double terminal_cross_track_weight = 100.0;
  double terminal_heading_weight = 8.0;
  double terminal_world_velocity_weight = 8.0;
  double terminal_yaw_rate_weight = 2.0;
};

// Wheel-rate-constrained dynamic LTV MPC. The prediction state contains pose,
// world-frame velocity, and yaw rate. At every solve, the stateful drift model
// is linearized around the reference trajectory and condensed into a QP.
class LtvMpc {
 public:
  explicit LtvMpc(MpcConfig config);

  Twist2d Solve(const RobotState& state, const Trajectory& reference, double reference_time);
  const std::vector<RobotState>& predicted_states() const { return predicted_states_; }
  double last_cost() const { return last_cost_; }

 private:
  using ControlSequence = std::vector<Twist2d>;
  using DynamicState = Eigen::Matrix<double, 6, 1>;

  ControlSequence InitialGuess(const RobotState& state);
  void ProjectControls(ControlSequence* controls, Twist2d previous_control) const;
  std::vector<RobotState> Rollout(const RobotState& state, const ControlSequence& controls) const;
  DynamicState ToDynamicState(const RobotState& state) const;
  RobotState ToRobotState(const DynamicState& state) const;
  DynamicState StepModel(const DynamicState& state, const Eigen::Vector2d& control) const;
  DynamicState StateDifference(const DynamicState& left, const DynamicState& right) const;
  void LinearizeDynamics(const DynamicState& state, const Eigen::Vector2d& control,
                         Eigen::Matrix<double, 6, 6>* state_jacobian,
                         Eigen::Matrix<double, 6, 2>* control_jacobian) const;
  void BuildQuadratic(const RobotState& state, const Trajectory& reference,
                      double reference_time, Eigen::MatrixXd* hessian,
                      Eigen::VectorXd* gradient) const;
  double QuadraticCost(const ControlSequence& controls, const Eigen::MatrixXd& hessian,
                       const Eigen::VectorXd& gradient) const;

  MpcConfig config_;
  ControlSequence previous_solution_;
  Twist2d last_command_;
  std::vector<RobotState> predicted_states_;
  double last_cost_ = 0.0;
};

// Kept as a source-compatible name for callers of the original example. The
// implementation is now LTV MPC, not the old nonlinear finite-difference MPC.
using NonlinearMpc = LtvMpc;

}  // namespace differential_drive_mpc

#endif  // DIFFERENTIAL_DRIVE_MPC_NONLINEAR_MPC_H_
