
set(SOURCES
    cloth.cpp)

add_library(phys_lib ${SOURCES})
target_link_libraries(phys_lib INTERFACE solvant)

#add_subdirectory(tests)
