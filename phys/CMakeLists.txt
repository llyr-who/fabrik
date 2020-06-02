
set(SOURCES
    cloth.cpp)

add_library(object_lib ${SOURCES})
target_link_libraries(object_libraries ${OPENGL_LIBRARIES} ${GLFW_LIBRARIES})

add_subdirectory(tests)
