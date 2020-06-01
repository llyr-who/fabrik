cmake_minimum_required(VERSION 3.10)
project(fabrik VERSION 1.0 LANGUAGES CXX)

find_package(OpenGL REQUIRED)
find_package(X11 REQUIRED)

set(SOURCES 
 antmath.cpp
 cloth.cpp
 GeometryGen.cpp
 GLnixAPP.cpp
 GLnixDEMO2.cpp
 MathHelper.cpp
 Timer.cpp )

add_library(fabrik_lib ${SOURCES})
add_executable(fabrik_run GLnixDEMO2.cpp)
target_link_libraries(fabrik_run PRIVATE fabrik_lib 
                                 ${OPENGL_LIBRARIES} 
                                 ${X11_LIBRARIES}
                                 ${X11_Xrender_LIB})

