#
# set options
#
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

if(NOT CMAKE_CXX_COMILER)
  set(CMAKE_CXX_COMILER mpicxx)
endif()

if(NOT BOOST_DIR)
  message(FATAL_ERROR "required -DBOOST_DIR=<path of boost v1.75> option")
endif()

option(OPTIMIZE_FOR_NATIVE "Build with -march=native" ON)
if(OPTIMIZE_FOR_NATIVE AND ${CMAKE_BUILD_TYPE} STREQUAL Release)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  message("add -merge=native to CMAKE_CXX_FLAGS")
endif()

option(SHARED_MEMORY_ONLY "Compile only shared-memory version" OFF)
