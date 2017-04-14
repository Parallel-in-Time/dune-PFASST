if(NOT dune-PFASST_FOUND)
# Whether this module is installed or not
set(dune-PFASST_INSTALLED OFF)

# Settings specific to the module

# Package initialization
# Set prefix to source dir
set(PACKAGE_PREFIX_DIR /home/zam/ruth/software_engineering/dune/dune-PFASST)
macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

#report other information
set_and_check(dune-PFASST_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(dune-PFASST_INCLUDE_DIRS "/home/zam/ruth/software_engineering/dune/dune-PFASST")
set(dune-PFASST_CXX_FLAGS " -std=c++14 ")
set(dune-PFASST_CXX_FLAGS_DEBUG "-g")
set(dune-PFASST_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(dune-PFASST_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(dune-PFASST_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(dune-PFASST_DEPENDS "dune-common;dune-geometry;dune-grid;dune-localfunctions;dune-istl;dune-typetree;dune-functions;dune-matrix-vector;dune-fufem")
set(dune-PFASST_SUGGESTS "")
set(dune-PFASST_MODULE_PATH "/home/zam/ruth/software_engineering/dune/dune-PFASST/cmake/modules")
set(dune-PFASST_LIBRARIES "")

# Lines that are set by the CMake build system via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


#import the target
if(dune-PFASST_LIBRARIES)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/dune-PFASST-targets.cmake")
endif()
endif()
