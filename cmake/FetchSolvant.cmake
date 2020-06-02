include(ExternalProject)

ExternalProject_Add(
  solvant
  GIT_REPOSITORY https://github.com/llyr-who/solvant
  GIT_TAG origin/master 
  )

ExternalProject_Get_Property(solvant install_dir)

include_directories(${install_dir}/src/solvant)
