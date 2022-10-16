# - Find GMP
# Find the GMP includes and libraries
# This module defines
#  GMP_INCLUDE_DIR, where to find gmp.h
#  GMP_LIBRARY, libgmp.so for UNIX
#  GMP_FOUND, If false, do not try to use GMP.
#
# Available targets:
#  GMP::GMP

find_path(GMP_INCLUDE_DIR gmp.h
  PATHS
    ENV GMP_ROOT
    ENV GMP_INCLUDE_DIR
    ${GMP_ROOT}
    /usr
    /usr/local
    $ENV{HOME}/.local
  PATH_SUFFIXES
    include
  )

find_library(GMP_LIBRARY
  NAMES
    gmp
  PATHS
    ENV GMP_ROOT
    ENV GMP_LIB_DIR
    ${GMP_ROOT}
    /usr
    /usr/local
  PATH_SUFFIXES
    lib
  )

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
  REQUIRED_VARS
    GMP_INCLUDE_DIR
    GMP_LIBRARY
  )

if(GMP_FOUND AND NOT TARGET GMP::GMP)
  add_library(GMP::GMP UNKNOWN IMPORTED)
  set_target_properties(GMP::GMP PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES "C"
    IMPORTED_LOCATION "${GMP_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
    )
endif()
