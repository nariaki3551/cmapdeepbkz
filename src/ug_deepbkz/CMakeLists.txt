cmake_minimum_required(VERSION 3.1)
project(
  CMAP-DeepBKZ
  DESCRIPTION "Massively parallelized DeepBKZ based on CMAP-LAP"
  LANGUAGES CXX
  )


#
# synchronous_cmapdeepbkz -- openMP shared memory parallel version with synchronous
#
add_executable(
    synchronous_cmapdeepbkz
    synchronousCmapDeepBkz.cpp
    )
target_link_libraries(
    synchronous_cmapdeepbkz PRIVATE
    laptoolslib
    ${BASE_LIB}
    )
target_include_directories(
    synchronous_cmapdeepbkz PRIVATE
    ${CMAKE_SOURCE_DIR}/src/laptools
    )


#
# synchronous_cmapdeepbkz_identity -- openMP shared memory parallel version with synchronous using identity matrix
#
add_executable(
    synchronous_cmapdeepbkz_identity
    synchronousCmapDeepBkzIdentity.cpp
    )
target_link_libraries(
    synchronous_cmapdeepbkz_identity PRIVATE
    laptoolslib
    ${BASE_LIB}
    )
target_include_directories(
    synchronous_cmapdeepbkz_identity PRIVATE
    ${CMAKE_SOURCE_DIR}/src/laptools
    )


#
# fcmapdeepbkz -- shared memory parallel version
#
add_executable(
    fcmapdeepbkz
    fcmapdeepbkz.cpp
    cmapDeepBkzParaLoadCoordinator.cpp
    cmapDeepBkzParaSolver.cpp
    cmapDeepLll.cpp
    cmapDeepBkz.cpp
    )
target_link_libraries(
    fcmapdeepbkz PRIVATE
    fcmaplaplib
    laptoolslib
    ${BASE_LIB}
    )
target_include_directories(
    fcmapdeepbkz PRIVATE
    ${CMAKE_SOURCE_DIR}/src/ug_cmaplap
    ${CMAKE_SOURCE_DIR}/src/laptools
    )


#
# paracmapdeepbkz -- distributed memory parallel version
#
# set(CMAKE_CXX_COMPILER mpicxx)
# set(CMAKE_CXX_COMPILER mpiicpc)
if( SHARED_MEMORY_ONLY )
  return()
endif()

find_package( MPI )
add_executable(
    paracmapdeepbkz
    paracmapdeepbkz.cpp
    cmapDeepBkzParaLoadCoordinator.cpp
    cmapDeepBkzParaSolver.cpp
    cmapDeepLll.cpp
    cmapDeepBkz.cpp
    )
target_link_libraries(
    paracmapdeepbkz PRIVATE
    paracmaplaplib
    laptoolslib
    ${BASE_LIB}
    )
target_include_directories(
    paracmapdeepbkz PRIVATE
    ${CMAKE_SOURCE_DIR}/src/ug_cmaplap
    ${CMAKE_SOURCE_DIR}/src/laptools
    )
if( MPI_COMPILE_FLAGS )
    set_target_properties(
        paracmapdeepbkz PROPERTIES
        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
        )
endif()
