cmake_minimum_required(VERSION 3.27.4)
project(time_series)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(time_series SHARED TimeSeries.cpp)
set_target_properties(time_series PROPERTIES PUBLIC_HEADER TimeSeries.h)#Used to define this as a a shared library (meant to be used by multiple programs)

install(TARGETS time_series LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib
        PUBLIC_HEADER DESTINATION ${CMAKE_SOURCE_DIR}/include)

