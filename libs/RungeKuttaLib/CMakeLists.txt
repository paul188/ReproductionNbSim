cmake_minimum_required(VERSION 3.27.4)
project(runge_kutta)

set(${Eigen3_DIR} ../../../../../../../Eigen_install)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(runge_kutta STATIC RungeKutta.cpp)
set_target_properties(runge_kutta PROPERTIES PUBLIC_HEADER RungeKutta.h)#Used to define this as a a shared library (meant to be used by multiple programs)

target_link_libraries(runge_kutta ../libs/TimeSeriesLib/build/Debug/time_series.lib
../libs/ParamsLib/build/Debug/params.lib)

install(TARGETS runge_kutta LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include/runge_kutta)

