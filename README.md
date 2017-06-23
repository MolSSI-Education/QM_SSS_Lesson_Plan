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



