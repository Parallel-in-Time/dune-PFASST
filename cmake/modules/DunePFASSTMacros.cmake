# File for module specific CMake tests.

#find_package(Eigen3 REQUIRED)
find_package(Eigen3)


if(EIGEN3_FOUND)
  dune_register_package_flags(INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
endif()

function(dune_add_pfasst_flags _targets)
  if(EIGEN3_FOUND)
    set_property(TARGET ${_targets_} APPEND PROPERTY INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
  endif()
endfunction()

#dune-PFASST needs an installation of PFASST++
#set installation path with cmake -Dpfasst_include=...
#SET(pfasst_include="/home/zam/ruth/software_engineering/PFASST/include" CACHE STRING "PFASST_DIR")#change default directory to /usr/include before push
#include_directories(${pfasst_include})
#MESSAGE( STATUS "include directory of PFASST++ library is set to " ${pfasst_include} )
