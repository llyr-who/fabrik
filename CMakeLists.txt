cmake_minimum_required(VERSION 3.10)
project(fabrik VERSION 1.0 LANGUAGES CXX)

find_package(OpenGL REQUIRED)
find_package(glfw 3.0 REQUIRED)

include(cmake/FetchSolvant.cmake)

#add_subdirectory(geo)
add_subdirectory(phys)
#add_subdirectory(utils)
