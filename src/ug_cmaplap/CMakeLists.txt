cmake_minimum_required(VERSION 3.1)
project(
    CMAP-LAP
    DESCRIPTION "Massively Parallel Framework for Lattice Problems"
    LANGUAGES CXX
    )

#
# setting for fcmaplaplib and paracmaplaplib
#
set(
    UG_BASE
    ../ug/paraLoadCoordinator.cpp
    ../ug/paraParamSet.cpp
    ../ug/paraSolver.cpp
    ../ug/gzstream.cpp
    ../ug/uggithash.cpp
    ../ug/paraSysTimer.cpp
    ../ug/githash.cpp
    )
set(
    CMAPLAP_BASE
    cmapLapParaBasis.cpp
    cmapLapParaBasisTh.cpp
    cmapLapParaCalculationStateTh.cpp
    cmapLapParaCheckpointWriter.cpp
    cmapLapParaInitiator.cpp
    cmapLapParaInstanceTh.cpp
    cmapLapParaLattice.cpp
    cmapLapParaLoadCoordinator.cpp
    cmapLapParaLoadCoordinatorTerminationState.cpp
    cmapLapParaPackedVectorTh.cpp
    cmapLapParaParamSet.cpp
    cmapLapParaShareDataPool.cpp
    cmapLapParaSolution.cpp
    cmapLapParaSolutionTh.cpp
    cmapLapParaSolver.cpp
    cmapLapParaSolverLocalComm.cpp
    cmapLapParaSolverPool.cpp
    cmapLapParaSolverStateTh.cpp
    cmapLapParaSolverTerminationState.cpp
    cmapLapParaSolverTerminationStateTh.cpp
    cmapLapParaTask.cpp
    cmapLapParaTaskTh.cpp
    )


#
# fuglib -- shared memory ug library
#
add_library(
    fuglib STATIC
    ${UG_BASE}
    ../ug/paraCommCPP11.cpp
    )
target_compile_definitions(
    fuglib PUBLIC
    -DUG_WITH_ZLIB
    -D_COMM_CPP11
    )
target_link_libraries(
    fuglib
    ${BASE_LIB}
    )

#
# parauglib -- distributed memory ug library
#
if( NOT SHARED_MEMORY_ONLY )
add_library(
    parauglib STATIC
    ${UG_BASE}
    ../ug/paraCommMpi.cpp
    ../ug/paraParamSetMpi.cpp
    ../ug/paraTimerMpi.cpp
    )
target_include_directories(
    parauglib PUBLIC
    ${MPI_INCLUDE_PATH}
    )
target_compile_definitions(
    parauglib PUBLIC
    -DUG_WITH_ZLIB
    -D_COMM_MPI_WORLD
    -DMPICH_IGNORE_CXX_SEEK
    )
target_link_libraries(
    parauglib
    ${BASE_LIB}
    )
endif()



#
# fcmaplaplib -- shared memory parallel version
#
add_library(
    fcmaplaplib STATIC
    ${CMAPLAP_BASE}
    cmapLapParaCommTh.cpp
    )
target_include_directories(
    fcmaplaplib PUBLIC
    ${CMAKE_SOURCE_DIR}/usr/include
    )
target_compile_definitions(
    fcmaplaplib PUBLIC
    -DUG_WITH_ZLIB
    -D_COMM_CPP11
    )
target_link_libraries(
    fcmaplaplib
    fuglib
    ${BASE_LIB}
    )


#
# paracmaplaplib -- distributed memory parallel version
#
# set(CMAKE_CXX_COMPILER mpicxx)
# set(CMAKE_CXX_COMPILER mpiicpc)
if( SHARED_MEMORY_ONLY )
  return()
endif()

find_package( MPI )

add_library(
    paracmaplaplib STATIC
    ${CMAPLAP_BASE}
    cmapLapParaCommMpi.cpp
    cmapLapParaInstanceMpi.cpp
    cmapLapParaSolutionMpi.cpp
    cmapLapParaTaskMpi.cpp
    cmapLapParaSolverStateMpi.cpp
    cmapLapParaCalculationStateMpi.cpp
    cmapLapParaSolverTerminationStateMpi.cpp
    cmapLapParaPackedVectorMpi.cpp
    cmapLapParaBasisMpi.cpp
    )
target_include_directories(
    paracmaplaplib PRIVATE
    ${CMAKE_SOURCE_DIR}/usr/include
    ${MPI_INCLUDE_PATH}
    )
target_compile_definitions(
    paracmaplaplib PUBLIC
    -DUG_WITH_ZLIB
    -D_COMM_MPI_WORLD
    -DMPICH_IGNORE_CXX_SEEK
    )
target_link_libraries(
    paracmaplaplib
    parauglib
    ${BASE_LIB}
    )

