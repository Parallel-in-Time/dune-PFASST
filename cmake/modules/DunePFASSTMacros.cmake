# File for module specific CMake tests.

find_package(Eigen3)


if(EIGEN3_FOUND)
  dune_register_package_flags(INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
endif()

function(dune_add_pfasst_flags _targets)
  if(EIGEN3_FOUND)
    set_property(TARGET ${_targets_} APPEND PROPERTY INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
  endif()
endfunction()
