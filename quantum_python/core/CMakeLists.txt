cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(core)

set(PYBIND11_CPP_STANDARD -std=c++11)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(core core.cc)

find_package(TargetLAPACK REQUIRED)
target_link_libraries(core PRIVATE tgt::lapack)
