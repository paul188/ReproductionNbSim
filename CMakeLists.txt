cmake_minimum_required(VERSION 3.27.4)

project(model)

add_executable(model ./src/cpp/main.cpp)

target_link_libraries(model ${CMAKE_SOURCE_DIR}/libs/TimeSeriesLib/build/Debug/time_series.lib
${CMAKE_SOURCE_DIR}/libs/RungeKuttaLib/build/Debug/runge_kutta.lib
${CMAKE_SOURCE_DIR}/libs/ParamsLib/build/Debug/params.lib)