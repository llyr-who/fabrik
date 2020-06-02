include(FetchContent)

FetchContent_Declare(
  solvant
  GIT_REPOSITORY https://github.com/llyr-who/solvant.git
  GIT_TAG        master
)

FetchContent_MakeAvailable(solvant)

message(STATUS ${solvant_SOURCE_DIR})
