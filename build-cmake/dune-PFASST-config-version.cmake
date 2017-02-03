set(PACKAGE_VERSION "1.0.0")

if("${PACKAGE_FIND_VERSION_MAJOR}" EQUAL "1" AND
     "${PACKAGE_FIND_VERSION_MINOR}" EQUAL "0")
  set (PACKAGE_VERSION_COMPATIBLE 1) # compatible with newer
  if ("${PACKAGE_FIND_VERSION}" VERSION_EQUAL "1.0.0")
    set(PACKAGE_VERSION_EXACT 1) #exact match for this version
  endif()
endif()
