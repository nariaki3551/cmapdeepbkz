cmake_minimum_required(VERSION 3.1)
project(
  cmaptest
  DESCRIPTION "Test configure of CMAP-LAP framework"
  LANGUAGES CXX
  )


#
# fcmaptest -- shared memory parallel version
#

add_executable(
    fcmaptest
    fcmaptest.cpp
    cmapTestParaSolver.cpp
    cmapTestDeepBkz.cpp
    cmapTestDeepLll.cpp
    cmapTestEnumeration.cpp
    cmapTestGaussSieve.cpp
    )
target_link_libraries(
    fcmaptest PRIVATE
    fcmaplaplib
    laptoolslib
    ${BASE_LIB}
    )
target_include_directories(
    fcmaptest PRIVATE
    ${CMAKE_SOURCE_DIR}/src/ug_cmaplap
    ${CMAKE_SOURCE_DIR}/src/laptools
    )


#
# paracmaptest -- distributed memory parallel version
#
# set(CMAKE_CXX_COMPILER mpicxx)
# set(CMAKE_CXX_COMPILER mpiicpc)
if( SHARED_MEMORY_ONLY )
  return()
endif()

find_package( MPI )
add_executable(
    paracmaptest
    paracmaptest.cpp
    cmapTestParaSolver.cpp
    cmapTestDeepBkz.cpp
    cmapTestDeepLll.cpp
    cmapTestEnumeration.cpp
    cmapTestGaussSieve.cpp
    )
target_link_libraries(
    paracmaptest PRIVATE
    paracmaplaplib
    laptoolslib
    ${BASE_LIB}
    )
target_include_directories(
    paracmaptest PRIVATE
    ${CMAKE_SOURCE_DIR}/src/ug_cmaplap
    ${CMAKE_SOURCE_DIR}/src/laptools
    )
if( MPI_COMPILE_FLAGS )
    set_target_properties(
        paracmaptest PROPERTIES
        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
        )
endif()

