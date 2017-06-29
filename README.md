<p align="center">

<!--Travis-->
<a href="https://travis-ci.org/dgasmith/QM_SSS_Lesson_Plan"><img src="https://travis-ci.org/dgasmith/QM_SSS_Lesson_Plan.svg?branch=master"></a>

<!--CodeCov-->
<a href="https://codecov.io/gh/dgasmith/QM_SSS_Lesson_Plan">
  <img src="https://codecov.io/gh/dgasmith/QM_SSS_Lesson_Plan/branch/master/graph/badge.svg" alt="Codecov" />
</a>

<!--License-->
<a href="https://opensource.org/licenses/LGPL-3.0"> <img src="https://img.shields.io/github/license/dgasmith/QM_SSS_Lesson_Plan.svg" /></a>

</p>

# QM_SSS_Lesson_Plan
This repository will provide a starting template for new Python-based QM projects.
Here we pick a set of tools that we are familiar with; however, there are many
tools extant that have similar functionality. 

The following tools will be used:
 - [GitHub](github.com) - Version Control
 - [Travis CI](https://travis-ci.org) - Continous Integration
 - [pytest](https://docs.pytest.org/en/latest/) - Unit and Regression Testing
 - [CodeCov](https://codecov.io) - Testing Coverage Analysis
 - [Read the Docs](https://readthedocs.org) - Documentation

## Testing

To run the test suite please first run `pip install -e .` in the base
repository folder. This will register this repository with your local Python so
that `import quantum_python` will work in any directory. Tests can then be run
with `py.test -v`. If `pytest` is not found please use `pip install pytest` to
aquire the module.



