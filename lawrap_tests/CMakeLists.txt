cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(lawrap_tests)

add_executable(test_c test.c)
add_executable(test_cxx test.cc)

find_package(LAWrap REQUIRED)
target_link_libraries(test_c PRIVATE LAWrap::lawrap)
target_link_libraries(test_cxx PRIVATE LAWrap::lawrap)

# Usage:
#  >>> ls
#  CMakeLists.txt   test.c    test.cc
#  >>> cmake -H. -Bbuild
#  -- The C compiler identification is GNU 5.2.0
#  -- The CXX compiler identification is GNU 5.2.0
#  ...
#  -- TargetLAPACKConfig detected.
#  -- Found TargetLAPACK: Found LAPACK MKL: /home/psilocaluser/miniconda3/envs/sss/lib/libmkl_rt.so;/home/psilocaluser/miniconda3/envs/sss/lib/libiomp5.so;
#  -- Configuring done
#  -- Generating done
#  ...
#  >>> make -j2
#  Scanning dependencies of target test_cxx
#  ...
#  [100%] Built target test_cxx
#
#  >>> ./test_cxx
#  info 0

# Info:
# The CMake install for LAWrap does LAPACK detection
# at `find_package(LAWrap)` time, so linking against it
# (`target_link_libraries(tgt LAWrap::lawrap)`) will pull in
# BLAS/LAPACK. The detection is in three layers, in decreasing order
# of precedence:
#
# (1) libraries and optionally headers passed directly as ;-separated list. E.g.,
#     `cmake -DLAPACK_LIBRARIES="/path/libmkl_rt.so;/path/libiomp5.so" -DLAPACK_INCLUDE_DIRS=/path/mkl.h`
#
#
# (2) finds a TargetLAPACKConfig.cmake file that specifies LAPACK/BLAS
#     targets, including libraries and headers. This is provided by Psi4 and
#     is how the SSS is set up to work. Influence its detection analogously to
#     https://github.com/loriab/lawrap/blob/cmake/cmake/LAWrapConfig.cmake.in#L38-L42.
#
# (3) wraps ordinary cmake-provided `find_package(LAPACK)`. This
#     works just fine for libraries, but it doesn't detect headers. I could
#     easily set up the TargetLAPACK to find a mkl.h, but finding the headers
#     _corresponding_ to the found libraries would be quite tricky.
