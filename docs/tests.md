# Tests
We made preconstructed test with fixed parameters to stress test the library.

## Running Tests

- `test_solver` (baseline sin(πx))
- `test_solver_highfreq` (sin(3πx))
- `test_solver_reaction` (strong reaction term)
- `test_solver_polynomial` (cubic-like profile)
- `test_solver_shifted_domain` (domain [0,2])
- `test_solver_overlap_steep` (stiffer diffusion/reaction balance)

Run them individually from the project root folder.
For example, to run `test_solver` write in the terminal:

```bash
./bin/test_solver
```
Some tests produce outputs in the `outputs` directory.
These are `.vtk` files, for more information on how to visualize them see [Visualization](paraview_visualization.md)

