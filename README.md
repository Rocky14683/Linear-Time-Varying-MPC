# Differential-Drive MPC

This C++20 project simulates and controls a differential-drive robot that follows a
time-scaled 2D Bezier trajectory. It is inspired by the organization and visualization
workflow of [Rocky14683/MCL](https://github.com/Rocky14683/MCL), but implements
trajectory tracking rather than localization.

## What is included

- `BezierTrajectory`: analytic position, velocity, acceleration, heading, linear/angular
  velocity, and linear/angular acceleration. It automatically lengthens the trajectory
  duration until every supplied chassis *and composed wheel* motion limit is met.
- `DifferentialDriveRobot`: converts linear/angular commands to left/right wheel speeds,
  simulates pose integration, and supports position-varying friction, turn-induced lateral
  drift, plus world-frame external acceleration and angular acceleration.
- `LtvMpc`: constrained, warm-started dynamic linear time-varying MPC. Its prediction state
  contains pose, world-frame velocity, and yaw rate, and it condenses the locally linearized
  drift model into a QP under wheel-rate and wheel-acceleration bounds.
- `RerunVisualizer`: logs the target curve, MPC-followed path, robot footprint and heading,
  MPC horizon, tracking errors, and velocity/wheel telemetry to Rerun.
- GoogleTest unit tests covering trajectory feasibility, drive kinematics/disturbances,
  MPC command constraints, and closed-loop cross-track correction.

## Build and run

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/run_mpc_simulation
```

The simulation opens Rerun and streams a live, time-indexed view. The target trajectory is
blue, with its visible dots sampled uniformly by arc length; the MPC-followed path is red, the
simulated robot is orange, and the current MPC prediction is green. Bezier control points are
yellow, joined by a purple control polygon.
The MPC-followed path is a linear-velocity heatmap: blue is stopped, cyan/green is moderate,
and yellow/red approaches the configured linear-velocity limit. The timeline is named
`simulation_time`. Each plot is logged as one multi-series `Scalars` batch: linear velocity,
angular velocity, wheel speed, or tracking error.

## Design notes

The reference heading is the tangent angle of the Bezier curve. Therefore it satisfies the
differential-drive nonholonomic constraint: the world-frame velocity is always aligned with
the robot heading. The generator first builds an arc-length lookup table, then applies the
minimum-jerk profile to distance \(d\), rather than to the Bezier control parameter. Thus
linear velocity is \(\dot{d}\) regardless of how quickly the Bezier parameter changes. It
rejects cusps, where the heading would be undefined, and uses analytic derivatives through
third order to constrain both linear and angular acceleration.

The MPC linearizes the same stateful drift dynamics used by the plant at each reference-horizon
sample. Its state is pose, world-frame velocity, and yaw rate, so it predicts lateral momentum
through a turn and can command corrective yaw before that momentum becomes cross-track error.
The physical wheel constraints are enforced in wheel coordinates:
`wheel_rate = (v +/- omega * track_width / 2) / wheel_radius`. This is a diamond-shaped coupled
constraint in `(v, omega)` space, so independently limiting `v` and `omega` cannot substitute
for it.

`DifferentialDriveConfig::lateral_drift_coefficient` controls turn-induced inertia in the plant.
The plant first integrates drive and external acceleration into `RobotState::world_velocity`, then
retains this fraction of the velocity component lateral to the new heading. `0.0` enforces ideal
no-slip kinematics; `1.0` keeps all lateral inertia. This gives drift memory across time steps
instead of an instantaneous lateral-velocity offset.

## Slip-based prediction model

The dynamic MPC and plant use the same normalized two-side contact model. Resolve world velocity
into the body frame as `u = heading · v_world` and `v = lateral · v_world`. The wheel-surface-speed
mismatches are

\[
s_L=(v_c-\tfrac b2\omega_c)-(u-\tfrac b2r),\qquad
s_R=(v_c+\tfrac b2\omega_c)-(u+\tfrac b2r).
\]

The longitudinal acceleration of each side is proportional to its slip and capped by traction.
The available lateral acceleration is reduced through a friction circle:

\[
a_{x,i}=\operatorname{sat}(\tfrac{k_x}{2}s_i,\,a_{\max}/2),\qquad
|a_{y,i}|\le\sqrt{(a_{\max}/2)^2-a_{x,i}^2}.
\]

Each side requests lateral acceleration `-k_y v / 2`. Thus forward or differential wheel command
consumes traction and reduces lateral correction; the MPC predicts and compensates that coupling.

With unit forward and lateral vectors `t=[cos(theta), sin(theta)]^T` and
`n=[-sin(theta), cos(theta)]^T`, the continuous state equations used before discretization are

\[
\dot p=V^W,\qquad
\dot V^W=t(a_{x,L}+a_{x,R})+n(a_{y,L}+a_{y,R}),\qquad
\dot\theta=r,
\]

\[
\dot r=\frac{k_r}{b}(s_R-s_L)-d_r r.
\]

After each Euler step, the model retains the configured fraction `c` of the component of
`V^W` lateral to the *new* heading. This is the accumulated-drift state: `c=0` is ideal
no-slip behavior and `c=1` preserves all lateral inertia. At every MPC update, this nonlinear
transition is finite-difference linearized about the latest predicted state/control trajectory
to produce the local LTV model used by the QP.

Position error is resolved in the reference-path frame: along-track error is parallel to the
reference tangent and cross-track error is normal to it. Cross-track error has a higher cost,
so the controller favors remaining on the target curve over matching its progress exactly.
Angular command magnitude is also penalized to discourage unnecessary turning.
