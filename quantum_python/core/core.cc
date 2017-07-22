#include <exception>
#include <stdlib.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <omp.h>

namespace py = pybind11;

void compute_PKJK(py::array_t<double> I, py::array_t<double> D, py::array_t<double> J, py::array_t<double> K) {

    // Grab the py::buffer info from each array
    py::buffer_info I_info = I.request();
    py::buffer_info D_info = D.request();
    py::buffer_info J_info = J.request();
    py::buffer_info K_info = K.request();

    // Do a few size checks
    size_t nbf = I_info.shape[0];
    size_t nbf2 = nbf * nbf;
    size_t nbf3 = nbf2 * nbf;

    if ((D_info.shape.size() != 2) || (nbf != D_info.shape[0]) || (nbf != D_info.shape[1])){
        throw std::length_error("The shape of the Density matrix does not match the ERI!");
    }

    if ((J_info.shape.size() != 2) || (nbf != J_info.shape[0]) || (nbf != J_info.shape[1])){
        throw std::length_error("The shape of the Coulomb matrix does not match the ERI!");
    }

    if ((K_info.shape.size() != 2) || (nbf != K_info.shape[0]) || (nbf != K_info.shape[1])){
        throw std::length_error("The shape of the Exchange matrix does not match the ERI!");
    }

    // Loop and build J and K
    double* I_ptr = static_cast<double*>(I_info.ptr);
    double* D_ptr = static_cast<double*>(D_info.ptr);
    double* J_ptr = static_cast<double*>(J_info.ptr);
    double* K_ptr = static_cast<double*>(K_info.ptr);

    for (size_t p = 0; p < nbf; p++){
        for (size_t q = 0; q < nbf; q++){
            for (size_t r = 0; r < nbf; r++){
                # pragma omp simd
                for (size_t s = 0; s < nbf; s++){
                    // printf("%zu %zu %zu %zu | %zu %zu\n", p, q, r, s, p * nbf + s, p * nbf3 + q * nbf2 + r * nbf + s);
                    J_ptr[p * nbf + q] += D_ptr[r * nbf + s] * I_ptr[p * nbf3 + q * nbf2 + r * nbf + s];
                    K_ptr[p * nbf + q] += D_ptr[r * nbf + s] * I_ptr[p * nbf3 + r * nbf2 + q * nbf + s];
                }
            }
        }
    }

}

PYBIND11_PLUGIN(core) {
    py::module m("core", "pybind11 core plugin");

    m.def("compute_PKJK", &compute_PKJK, "A function that can compute the PK J and K matrices.");

    return m.ptr();
}
