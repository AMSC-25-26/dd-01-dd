# Building the project

## Table of contents
1. [Requirements](#Requirements)
2. [How to compile](#How to compile)
3. [Execution](#Execution)

## Requirements

This project library compilation requires:

 - CMake ≥ 3.28
 - OpenMP
 - MPI (in the future)
 - Eigen 3 (in the future)

If you intend to compile the docs too please see [Doxygen documentation](documentation.md) for additional requirements.

## How to compile

Before compilation run cmake as follows while being in the project root directory:
```bash
cmake -S ./ -B build
```
To compile move into the `build` directory and run
```bash
make
```

Note that you may choose what to compile with the following options:
 * **all** compile everything, including tests (but not the docs)
 * **library** only compile the library files
 * **tests** compile all the tests (this will also compile the library)
 * **library_file_name** only compile that cpp file (note: do not include the extension)
 * **test_name** only compile that test
 * **docs** compile that docs, see [Doxygen documentation](documentation.md) for more.

Examples are given below:
```bash
make docs
make all
make tests
```

## Execution
The Library does not provide a main executable, it's meant for the user to create its own cpp 
file and use this library.
The test will be compiled to `./bin` please see [Testing](tests.md) to see how to use them.
