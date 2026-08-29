#include "differential_drive_mpc/rerun_visualizer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace differential_drive_mpc {
namespace rr = rerun;
namespace {

rr::Position2D ToRerunPoint(const Eigen::Vector2d &point) {
    return {static_cast<float>(point.x()), static_cast<float>(point.y())};
}

std::vector <rr::Position2D> ToRerunPoints(const std::vector <RobotState> &states) {
    std::vector <rr::Position2D> result;
    result.reserve(states.size());
    for (const RobotState &state: states) {
        result.push_back(ToRerunPoint(state.pose.position));
    }
    return result;
}

std::vector <rr::Position2D> ToRerunPoints(const std::vector <Eigen::Vector2d> &points) {
    std::vector <rr::Position2D> result;
    result.reserve(points.size());
    for (const Eigen::Vector2d &point: points) {
        result.push_back(ToRerunPoint(point));
    }
    return result;
}

std::vector <rr::Position2D> RobotFootprint(const Pose2d &pose) {
    constexpr float kHalfLength = 0.18F;
    constexpr float kHalfWidth = 0.12F;
    const Eigen::Vector2d forward(std::cos(pose.heading), std::sin(pose.heading));
    const Eigen::Vector2d lateral(-forward.y(), forward.x());
    return {ToRerunPoint(pose.position + kHalfLength * forward + kHalfWidth * lateral),
            ToRerunPoint(pose.position + kHalfLength * forward - kHalfWidth * lateral),
            ToRerunPoint(pose.position - kHalfLength * forward - kHalfWidth * lateral),
            ToRerunPoint(pose.position - kHalfLength * forward + kHalfWidth * lateral),
            ToRerunPoint(pose.position + kHalfLength * forward + kHalfWidth * lateral)};
}

}  // namespace

RerunVisualizer::RerunVisualizer(double linear_velocity_color_limit)
        : recording_("differential_drive_mpc_visualization_v3"),
          linear_velocity_color_limit_(std::max(linear_velocity_color_limit, 1e-6)) {
    ConfigurePlotStyles();
}

void RerunVisualizer::SpawnViewer() { recording_.spawn().exit_on_failure(); }

void RerunVisualizer::LogReference(const std::vector <Eigen::Vector2d> &reference) {
    target_trajectory_ = ToRerunPoints(reference);
}

void RerunVisualizer::LogControlPoints(const std::vector <Eigen::Vector2d> &control_points) {
    control_points_ = ToRerunPoints(control_points);
}

void RerunVisualizer::LogStep(double time, const TrajectoryPoint &target, const RobotState &state,
                              Twist2d command, WheelVelocities wheel_velocities,
                              const std::vector <RobotState> &prediction) {
    recording_.set_time_duration_secs("simulation_time", time);
    recording_.log("/world/target_trajectory/line", rr::LineStrips2D({target_trajectory_})
            .with_colors(rr::Color(55, 180, 255))
            .with_radii(0.040F));
    recording_.log(
            "/world/target_trajectory/waypoints",
            rr::Points2D(target_trajectory_).with_colors(rr::Color(55, 180, 255)).with_radii(0.007F));
    recording_.log(
            "/world/bezier_control_polygon",
            rr::LineStrips2D({control_points_}).with_colors(rr::Color(180, 120, 255)).with_radii(0.020F));
    recording_.log(
            "/world/bezier_control_points",
            rr::Points2D(control_points_).with_colors(rr::Color(255, 220, 70)).with_radii(0.02F));
    const rr::Position2D robot_position = ToRerunPoint(state.pose.position);
    followed_path_.push_back(robot_position);
    recording_.log(
            "/world/robot_position",
            rr::Points2D(followed_path_).with_colors(rr::Color(225, 0, 70)).with_radii(0.005F));
    const Eigen::Vector2d forward(std::cos(state.pose.heading), std::sin(state.pose.heading));
    const rr::Position2D heading_tip = ToRerunPoint(state.pose.position + 0.30 * forward);

    recording_.log(
            "/world/robot/position",
            rr::Points2D({robot_position}).with_colors(rr::Color(255, 115, 70)).with_radii(0.10F));
    recording_.log("/world/robot/footprint", rr::LineStrips2D({RobotFootprint(state.pose)})
            .with_colors(rr::Color(255, 115, 70))
            .with_radii(0.015F));
    const std::vector <std::vector<rr::Position2D>> heading_line{{robot_position, heading_tip}};
    recording_.log(
            "/world/robot/heading",
            rr::LineStrips2D(heading_line).with_colors(rr::Color(255, 220, 70)).with_radii(0.025F));
    recording_.log("/world/mpc_prediction", rr::LineStrips2D({ToRerunPoints(prediction)})
            .with_colors(rr::Color(70, 255, 140))
            .with_radii(0.020F));

    const Eigen::Vector2d position_error = state.pose.position - target.pose.position;
    const double heading_error = WrapAngle(state.pose.heading - target.pose.heading);
    recording_.log("/linear_velocity",
                   rr::Scalars({command.linear, state.linear_velocity, target.twist.linear}));
    recording_.log("/angular_velocity",
                   rr::Scalars({command.angular, state.angular_velocity, target.twist.angular}));
    recording_.log("/wheel_velocity", rr::Scalars({wheel_velocities.left, wheel_velocities.right}));
    recording_.log("/tracking_error", rr::Scalars({position_error.norm(), heading_error,
                                                   position_error.x(), position_error.y()}));
}

void RerunVisualizer::ConfigurePlotStyles() {
    const auto log_style = [this](const char *path, std::vector <rr::Rgba32> colors,
                                  std::vector <rr::components::Name> names) {
        recording_.log_static(path, rr::SeriesLines().with_colors(colors).with_names(names));
    };
    log_style("/linear_velocity",
              {rr::Rgba32{55, 180, 255}, rr::Rgba32{255, 140, 60}, rr::Rgba32{70, 255, 140}},
              {"command", "actual", "reference"});
    log_style("/angular_velocity",
              {rr::Rgba32{55, 180, 255}, rr::Rgba32{255, 140, 60}, rr::Rgba32{70, 255, 140}},
              {"command", "actual", "reference"});
    log_style("/wheel_velocity", {rr::Rgba32{120, 210, 90}, rr::Rgba32{205, 120, 255}},
              {"left wheel", "right wheel"});
    log_style("/tracking_error",
              {rr::Rgba32{255, 80, 80}, rr::Rgba32{255, 220, 70}, rr::Rgba32{80, 210, 255},
               rr::Rgba32{255, 120, 200}},
              {"position error", "heading error", "x error", "y error"});
}

}  // namespace differential_drive_mpc
