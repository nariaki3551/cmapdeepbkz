# - Find MPFR
# Find the MPFR includes and libraries
# This module defines
#  MPFR_INCLUDE_DIR, where to find mpfr.h
#  MPFR_LIBRARIES, the libraries needed to use MPFR.
#  MPFR_FOUND, If false, do not try to use MPFR.
#
# Available targets:
#  MPFR::MPFR

find_path(MPFR_INCLUDE_DIR mpfr.h
  PATHS
    ENV MPFR_ROOT
    ENV MPFR_INCLUDE_DIR
    ${MPFR_ROOT}
    /usr
    /usr/local
    $ENV{HOME}/.local
  PATH_SUFFIXES
    include
  )

find_library(MPFR_LIBRARY
  NAMES
    mpfr
  PATHS
    ENV MPFR_ROOT
    ENV MPFR_LIB_DIR
    ${MPFR_ROOT}
    /usr
    /usr/local
  PATH_SUFFIXES
    lib
  )

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR
  REQUIRED_VARS
    MPFR_INCLUDE_DIR
    MPFR_LIBRARY
  )

if (NOT TARGET GMP::GMP)
  find_package(GMP REQUIRED)
endif()

if(MPFR_FOUND AND NOT TARGET MPFR::MPFR)
  add_library(MPFR::MPFR UNKNOWN IMPORTED)
  set_target_properties(MPFR::MPFR PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES "C"
    IMPORTED_LOCATION "${MPFR_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES GMP::GMP
    )
endif()

