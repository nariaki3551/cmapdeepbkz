cmake_minimum_required(VERSION 3.18)
project(CMAP-LAP CXX)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

include(cmake/options.cmake)
include(cmake/base.cmake)
include(cmake/packages.cmake)

add_subdirectory(src/laptools)
add_subdirectory(src/ug_cmaplap)
add_subdirectory(src/ug_cmaptest)
add_subdirectory(src/ug_deepbkz)
