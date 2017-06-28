cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(core)

find_package(pybind11)

pybind11_add_module(core core.cc)
