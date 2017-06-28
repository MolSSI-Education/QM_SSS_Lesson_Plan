#! /usr/bin/env python

"""
This is a setup script for install a python project with a CMake dependency. 

For simplicity this will use Psi4's Cache to make detection extremely simple
"""

import setuptools
import distutils
import distutils.spawn as ds
import os
import subprocess as sp

try:
    import psi4
except ImportError:
    raise ImportError("Cannot find Psi4. Please make sure Psi4 is installed so"
                      "that this setup script can detect Psi4's CMake dependency trees.\n"
                      "Please run the following conda command: 'conda install psi4 -c psi4'")

class cmake_build(distutils.cmd.Command):

    description = 'Build the nested CMake project'
    user_options = [
#      # The format is (long option, short option, description).
#        ('pylint-rcfile=', None, 'path to Pylint config file'),
    ]

    def initialize_options(self):
      """Set default values for options."""
      # Each user option must be listed here with their default value.
  
    def finalize_options(self):
      """Post-process options."""

    def run(self):

        # Find build directory (in-place)
        abspath = os.path.abspath(os.path.dirname(__file__))
        build_path = os.path.join(abspath, "quantum_python", "core")
        os.chdir(build_path)

        # Capture cmake command 
        print("Aquiring CMake cache...")
        output = sp.check_output(["psi4", "--plugin-compile"]).decode("UTF-8")
        if "cmake -C" not in output:
            raise Exception("Psi4 Cache Error. Output as follows:\n" + output.decode("UTF-8"))

        # Run CMake command
        print("Building CMake structures...")
        output = sp.check_output(output.strip().split()).decode("UTF-8")
        if "Build files have been" not in output:
            raise Exception("CMake error. Output as follows:\n" + output)

        # Run install
        print("Compiling...")
        output = sp.check_output(["make", "-j2"]).decode("UTF-8")
        if "[100%]" not in output:
            raise Exception("Build error. Output as follows:\n" + output)

if __name__ == "__main__":
    setuptools.setup(
        name='quantum_python',
        version="0.2.1",
        description='A small quantum program written in Python with a Psi4 backend',
        author='Daniel Smith',
        author_email='dgasmith@vt.edu',
        url="https://github.com/dgasmith/quantum_python",
        license='BSD-3C',
        packages=setuptools.find_packages(),
        install_requires=[
            'numpy>=1.7',
            'pytest>=3.0',
            'pytest-cov',
        ],
        extras_require={
            'docs': [
                'sphinx==1.2.3',  # autodoc was broken in 1.3.1
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
        },

        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
        ],
        zip_safe=True,
    cmdclass={
        'cmake': cmake_build,
    },
    )
