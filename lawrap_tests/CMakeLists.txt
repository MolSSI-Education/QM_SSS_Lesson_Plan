cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(lawrap_tests)

add_executable(test_c test.c)
add_executable(test_cxx test.cc)

find_package(LAWrap REQUIRED)
target_link_libraries(test_c PRIVATE LAWrap::lawrap)
target_link_libraries(test_cxx PRIVATE LAWrap::lawrap)

