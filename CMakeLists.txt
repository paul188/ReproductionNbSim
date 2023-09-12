cmake_minimum_required(VERSION 3.27.4)

project(model)

add_executable(model ./src/cpp/main.cpp)

target_link_libraries(model C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/libs/TimeSeriesLib/build/Debug/time_series.lib
C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/libs/RungeKuttaLib/build/Debug/runge_kutta.lib)
