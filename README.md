# Domain Decomposition (1D Overlapping Schwarz)

This project implements a 1D overlapping Schwarz solver for linear diffusion–reaction problems using a tridiagonal (Thomas) factorization on each subdomain. Multiple validation cases and utilities are provided to exercise the solver and visualize results.

## Features
- PDESolver and SubdomainSolver classes specialized for 1D line domains.
- Overlapping Schwarz iterations with configurable subdomain count and overlap.
- Multiple analytic-test executables targeting different right-hand sides and boundary conditions.
- Optional VTK output for visualization (ParaView script in `src/visualization_pipeline.py`).

## Build
Requirements: CMake ≥ 3.28, C++17, OpenMP, Eigen3.

After cloning this repository, run the following
```bash
cd build
cmake ..
make
```

## Run Tests

- `test_solver` (baseline sin(πx))
- `test_solver_highfreq` (sin(3πx))
- `test_solver_reaction` (strong reaction term)
- `test_solver_polynomial` (cubic-like profile)
- `test_solver_shifted_domain` (domain [0,2])
- `test_solver_overlap_steep` (stiffer diffusion/reaction balance)

Run them individually from the `build` folder (current known limitation!). For example, to run `test_solver`

```bash
../bin/test_solver
```

Note that running this will produce an output file inside the `outputs` folder

## Results Visualization

Open [ParaView](https://www.paraview.org/) to visualize the VTK output files generated in the `outputs` folder. Use the provided `src/visualization_pipeline.py` script to visualize 1D solutions produced by the solver.

## Usage
To solve your own 1D problem, supply three small structs and a grid size:

1) `PDEParams` — coefficients and forcing
	- `mu` (diffusion), `c` (reaction)
	- `omega = {a, b}` domain limits
	- `dirichlet = {u_a, u_b}` boundary values at `a` and `b`
	- `f(x)` right-hand side (lambda or functor)

2) `SchwarzParams` — domain decomposition
	- `N` number of subdomains
	- `delta` overlap width (physical length)

3) `SolverParams` — iteration controls
	- `max_iter` maximum Schwarz iterations
	- `eps` stopping tolerance on the sup-norm difference between iterates

4) Grid spacing `h = (omega.b - omega.a) / (N_nodes - 1)`; `N_nodes` is your chosen number of grid points.

Typical solve sequence (as in test/test_solver.cpp):
```cpp
PDEParams pde{mu, c, f, omega, dirichlet};
SchwarzParams sp{N_subdomains, delta};
SolverParams sol{eps, max_iter};
Real h = (omega.b - omega.a) / (N_nodes - 1);

DiscreteSolver<Line> solver(pde, sp, &sol, h);
solver.solve();
auto u = solver.get_solution();

// Optional: write VTK for visualization
solver.print_to_file();
```

Notes:
- `N_subdomains * subdomain_width` should match the domain length; choose `delta` so overlaps are a few grid cells wide.
- `f` is evaluated pointwise on the grid the solver builds from `omega` and `h`.
- Check `solver.status` for convergence; `status.iter` reports iterations.
- Use `OMP_NUM_THREADS` to control OpenMP parallelism.
- VTK files are written to `outputs/`; view with `src/visualization_pipeline.py` in ParaView.

To use the solver for a custom problem, 

## Authors

| Name                   | PC       | GitHub Name       |
|------------------------|----------|-------------------|
| Samuele Allegranza     | 10766615 | samueleallegranza |
| Giulio Enzo Donninelli | 10823453 | gdonninelli       |
| Valeriia Potrebina     | 11114749 | ValeryPotrebina   |
| Alessia Rigoni         | 10859832 | alessiarigoni     |
| Vale Turco             | 10809855 | PurpleVale        |
