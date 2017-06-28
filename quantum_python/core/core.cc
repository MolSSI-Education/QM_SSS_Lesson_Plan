#include <stdlib.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}

PYBIND11_PLUGIN(core) {
    py::module m("core", "pybind11 core plugin");

    m.def("add", &add, "A function which adds two numbers");

    return m.ptr();
}
