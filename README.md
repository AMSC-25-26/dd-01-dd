Project Overview {#mainpage}
=========
This project implements a 1D overlapping Schwarz solver library for linear diffusion–reaction problems
using a tridiagonal (Thomas) factorization on each subdomain.
Multiple test cases and utilities are provided to try out the solver and visualize results.

The project is built with the vision of expansion in the future for the inclusion of
2D domains and cluster execution using MPI in conjunction with OpenMP.

# Table of contents
1. [Objective of the Projects](#Objective-of-the-Projects)
2. [How to read the docs](#How-to-read-the-docs)
3. [Project structure](#Project-structure)
4. [Class Structure](#Class-Structure)
5. [Example of Usage](#Example-of-Usage)
6. [Authors](#Authors)


# Objective of the Project
The objective is to create a solver for a 1 dimensional Elliptic PDE 
with constant diffusion and reaction coefficients.
We expect to expand it in the future to non-constant coefficients.

# How to read the docs
This page only contains overall documentation for the project, for in depth documentation
see the appropriate pages: 
- [Compilation and execution](docs/building.md)
- [Doxygen documentation](docs/documentation.md)
- [Testing](docs/tests.md)
- [Thomas Factorization](docs/thomas_algorithm.md)
- [Visualization](docs/paraview_visualization.md)

# Project structure
The project structure is as follows:
```
01-DD/
  ├── bin/
  ├── build/
  ├── docs/
  │   ├── Doxygen/
  │   └── pdfs/
  ├── include/
  ├── lib/
  ├── outputs/
  ├── src/
  └── test/
```
* `bin` will contain the executables (for example tests) created by cmake.
* `build` will be used by cmake, please see [Compilation and execution](docs/building.md) for more information.
* `docs/Doxygen` will contain the doxygen generated documentation, please see [Doxygen documentation](docs/documentation.md) for more information.
* `docs/pdfs` contains useful pds references pertaining DD methods.
* `include` contains all the`.hpp` files used by the library.
* `libs` will contain the compiled `.a` library files.
* `output` is a commodity directory used by test to save `.vtk` files.
* `src` contains the actual body of the library and a useful python executable used for visualization, see [Visualization](docs/paraview_visualization.md) for more information.
* `test` contains premade test to stress the library and ensure the correct behaviours of all its components.

# Class Structure
This project uses inheritance as a way to avoid repeated code, not much room was left for polymorphism.

All the classes are templates based on the `enum Dimension`
this is because we intend to extend the library to 2D in the future

Every class in the project inherits from `Types<dim>` which contains the definition for essential types.
For example in 1D the type `Domain` is simply a pair of 2 real numbers.

* The `PDESolver<dim>` class is not meant to be used, it's contains useful functions that are used by both
`SubdomainSolver<dim>` and `DiscreteSolver<dim>`.
It also contains attributes that must be known to both classes to correctly compute the solution.

* The `DiscreteSolver<dim>` class construct the approximation of the PDE and handles parallelism by dividing work and different subdomains.
It is the class to be called to actually obtain a solution, as it is in this class that we advance computation
and check for the correctness of the solution (based on the difference between iterates).

* The `SubdomainSolver<dim>` handles the construction of the algebraic system and finds the solution on the assigned subdomain.

* The `FactorizedTridiag<dim>` class contains the tridiagonal matrix and solves the system $Au = f$.
Before solving the system the matrix is factorised ($A = L*U$) we then solve 2 system in sequence:
$$\begin{cases} 
Ly=b \\
Ux = y
\end{cases}$$
We call `ftd.solve(b)` to solve the system on the r.h.s. `b`, please see [Thomas Algorithm](docs/thomas_algorithm.md) for more general info on the algorith used to solve the system.

For more information please read the compiled documentation.

# Example of Usage
To solve your own 1D problem, you need to create the following three small structs and select a grid size:

1) `PDEParams` — information about the real PDE problem
	- `mu` (diffusion coefficient), `c` (reaction coefficient)
	- `omega = {a, b}` domain limits
	- `dirichlet = {u_a, u_b}` boundary values at `a` and `b`
	- `f(x)` right-hand side (lambda or functor)

2) `SchwarzParams` — information about the domain decomposition parameters
	- `N` number of subdomains
	- `delta` overlap width (physical length)

3) `SolverParams` — iteration controls
	- `max_iter` maximum number of Schwarz iterations
	- `eps` stopping tolerance on the sup-norm difference between iterates

4) Grid spacing `h = (omega.b - omega.a) / (N_nodes - 1)`; `N_nodes` is your chosen number of grid points.

Typical usage to compute the solution is (as in test/test_solver.cpp):
```cpp
PDEParams pde{mu, c, f, omega, dirichlet};
SchwarzParams sp{N_subdomains, delta};
SolverParams sol{eps, max_iter};
Real h = (omega.b - omega.a) / (N_nodes - 1);

DiscreteSolver<Line> solver(pde, sp, sol, h);
solver.solve();
// check solver status (optional)
auto u = solver.get_solution();

// Optional: write VTK for visualization
solver.print_to_file(file);
```

Notes:
- `N_subdomains * subdomain_width` should match the domain length; choose `delta` so overlaps are a few grid cells wide.
- `f` is evaluated pointwise on the grid the solver builds from `omega` and `h`.
- Check `solver.status` for convergence; `solver.iter` reports iterations.
- Use `OMP_NUM_THREADS` environment variable to control OpenMP parallelism.
- VTK files are written to `outputs/`; view with `src/visualization_pipeline.py` in ParaView.

# Authors

| Name                   | PC       | GitHub Name       |
|------------------------|----------|-------------------|
| Samuele Allegranza     | 10766615 | samueleallegranza |
| Giulio Enzo Donninelli | 10823453 | gdonninelli       |
| Valeriia Potrebina     | 11114749 | ValeryPotrebina   |
| Alessia Rigoni         | 10859832 | alessiarigoni     |
| Vale Turco             | 10809855 | PurpleVale        |
