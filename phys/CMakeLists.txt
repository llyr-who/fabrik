set(HEADERS
    cloth.h)

set(SOURCES
    cloth.cpp)

add_library(phys_lib ${HEADERS} ${SOURCES})
target_link_libraries(phys_lib PUBLIC ${solvant_LIRARIES})
target_include_directories(phys_lib PUBLIC ${solvant_SOURCE_DIR})
#add_subdirectory(tests)
