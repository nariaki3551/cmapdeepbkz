cmake_minimum_required(VERSION 3.1)
project(
  laptools
  DESCRIPTION "Tools for lattice problems"
  LANGUAGES CXX
  )


add_library(
    laptoolslib
    Config.cpp
    Timer.cpp
    VectorElementPool.cpp
    Lattice.cpp
    Reduction.cpp
    DeepLll.cpp
    Enumeration.cpp
    DeepBkz.cpp
    KleinSampler.cpp
    GaussSieve.cpp
    SubEnumeration.cpp
    # RecursiveDeepBkz.cpp
    )
target_link_libraries(
    laptoolslib PRIVATE
    ${BASE_LIB}
    )

add_executable(
    seqcmaplap
    seqcmaplap.cpp
    )
target_link_libraries(
    seqcmaplap PRIVATE
    laptoolslib
    ${BASE_LIB}
    )


