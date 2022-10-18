#
# find third-party libraries
#
find_package( OpenMP REQUIRED )
find_package( ZLIB REQUIRED   )
find_library( NTL_LIB ntl           HINTS ${CMAKE_SOURCE_DIR}/usr/lib REQUIRED )
find_library( GMP_LIB gmp           HINTS ${CMAKE_SOURCE_DIR}/usr/lib REQUIRED )
find_library( PTHREAD_LIB pthread   HINTS ${CMAKE_SOURCE_DIR}/usr/lib REQUIRED )
find_package( BLAS )
find_package( MPFR REQUIRED )
find_package( GSL REQUIRED )

set(
    BASE_LIB
    OpenMP::OpenMP_CXX
    ZLIB::ZLIB
    ${NTL_LIB}
    ${GMP_LIB}
    ${PTHREAD_LIB}
    ${MPFR_LIBRARY}
    ${GSL_LIBRARY}
    )

if( BLAS_FOUND )
  add_definitions(
    -DEIGEN_USE_BLAS
    )
  message("add definitions -DEIGEN_USE_BLAS")

  list(
    APPEND
    BASE_LIB
    BLAS::BLAS
    )
endif()

