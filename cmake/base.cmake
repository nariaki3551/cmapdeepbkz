include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/usr/include
    )

if( CMAKE_BUILD_TYPE STREQUAL "Release")
    add_definitions(
        -DNDEBUG -DEIGEN_NO_DEBUG
        )
endif()


add_compile_options(
    -Wall
    -pedantic     # Warn of syntax that does not exist in ISO C/C++
    # -Wextra
    # -Wcast-align  # Warn of casts with larger alignment lengths
    # -Wcast-qual   # Warn for type modifiers in case of outliers
    # -Wctor-dtor-privacy
    # -Wformat=2
    -Winit-self   # Warn when an undefined variable initializes itself, such as int i = i
    -Wlogical-op  # Warn when logical operators are used where bitwise operators are appropriate
    # -Wmissing-declarations
    # -Wmissing-include-dirs
    -Wnoexcept
    # -Wold-style-cast
    # -Woverloaded-virtual
    -Wredundant-decls
    -Wshadow
    # -Wsign-conversion
    -Wsign-promo
    -Wstrict-null-sentinel
    -Wswitch-default
    -Wundef
    )
if( CMAKE_BUILD_TYPE STREQUAL "Release" )
    add_compile_options(
        -O3
        )
elseif( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    add_compile_options(
        -g
        )
endif()

