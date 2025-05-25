cmake_minimum_required(VERSION 3.10)

# Build gtest from sources in /usr/src/googletest
add_subdirectory(/usr/src/googletest googletest-build EXCLUDE_FROM_ALL)

# Compile Fortran code
add_library(t96_fortran STATIC T96.f)

# The C++ port library
add_library(t96_cpp_lib SHARED ../src/t96.cc)

# gtest test
add_executable(t96_compare testfortran_gtest.cpp)

target_include_directories(t96_compare PRIVATE ../src)

target_link_libraries(t96_compare
    PRIVATE t96_cpp_lib t96_fortran gtest_main)

add_test(NAME t96_compare_tests COMMAND t96_compare)
