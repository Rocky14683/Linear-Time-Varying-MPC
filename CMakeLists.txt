cmake_minimum_required(VERSION 3.25)

project(differential_drive_mpc LANGUAGES CXX)

option(BUILD_VISUALIZATION "Build the Rerun simulation executable." ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED NO_MODULE)

add_library(differential_drive_mpc
  src/bezier_trajectory.cc
  src/differential_drive_robot.cc
  src/nonlinear_mpc.cc
)

target_include_directories(differential_drive_mpc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(differential_drive_mpc PUBLIC Eigen3::Eigen)

target_compile_options(differential_drive_mpc PRIVATE
  -Wall
  -Wextra
  -Wpedantic
  -Werror
)

if(BUILD_VISUALIZATION)
  include(FetchContent)
  FetchContent_Declare(
    rerun_sdk
    URL https://github.com/rerun-io/rerun/releases/download/0.36.3/rerun_cpp_sdk.zip
  )
  FetchContent_MakeAvailable(rerun_sdk)

  add_executable(run_mpc_simulation
    app/run_mpc_simulation.cc
    src/rerun_visualizer.cc
  )
  target_link_libraries(run_mpc_simulation
    PRIVATE differential_drive_mpc rerun_sdk
  )
  target_compile_options(run_mpc_simulation PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()

include(CTest)
if(BUILD_TESTING)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  add_executable(differential_drive_mpc_test
    test/bezier_trajectory_test.cc
    test/differential_drive_robot_test.cc
    test/nonlinear_mpc_test.cc
  )
  target_link_libraries(differential_drive_mpc_test
    PRIVATE differential_drive_mpc GTest::gtest_main
  )
  target_compile_options(differential_drive_mpc_test PRIVATE -Wall -Wextra -Wpedantic -Werror)

  include(GoogleTest)
  gtest_discover_tests(differential_drive_mpc_test)
endif()
