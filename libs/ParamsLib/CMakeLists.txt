cmake_minimum_required(VERSION 3.27.4)
project(params)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(params SHARED params.cpp)
set_target_properties(params PROPERTIES PUBLIC_HEADER params.h)#Used to define this as a a shared library (meant to be used by multiple programs)

install(TARGETS params LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib
        PUBLIC_HEADER DESTINATION ${CMAKE_SOURCE_DIR}/include)

