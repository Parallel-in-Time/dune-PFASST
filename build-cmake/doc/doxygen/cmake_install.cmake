# Install script for directory: /home/zam/ruth/software_engineering/dune/dune-PFASST/doc/doxygen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/zam/ruth/software_engineering/dune_installation")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  execute_process(COMMAND /usr/bin/cmake --build /home/zam/ruth/software_engineering/dune/dune-PFASST/build-cmake --target doxygen_dune-PFASST
        WORKING_DIRECTORY /home/zam/ruth/software_engineering/dune/dune-PFASST/build-cmake/doc/doxygen)
      file(GLOB doxygenfiles
        GLOB /home/zam/ruth/software_engineering/dune/dune-PFASST/build-cmake/doc/doxygen/html/*.html
        /home/zam/ruth/software_engineering/dune/dune-PFASST/build-cmake/doc/doxygen/html/*.png
        /home/zam/ruth/software_engineering/dune/dune-PFASST/build-cmake/doc/doxygen/html/*.css
        /home/zam/ruth/software_engineering/dune/dune-PFASST/build-cmake/doc/doxygen/html/*.gif)
      set(doxygenfiles "${doxygenfiles}")
      foreach(_file ${doxygenfiles})
         get_filename_component(_basename ${_file} NAME)
         LIST(APPEND CMAKE_INSTALL_MANIFEST_FILES /home/zam/ruth/software_engineering/dune_installation/share/doc/dune-PFASST/doxygen/${_basename})
       endforeach()
       file(INSTALL ${doxygenfiles} DESTINATION /home/zam/ruth/software_engineering/dune_installation/share/doc/dune-PFASST/doxygen)
       message(STATUS "Installed doxygen into /home/zam/ruth/software_engineering/dune_installation/share/doc/dune-PFASST/doxygen")
endif()

