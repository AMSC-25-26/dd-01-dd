#include <types.hpp>
#include <solver.hpp>

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <filesystem>
#include <iomanip>
#include <omp.h>

using Real = Types<Line>::Real;
using Index = Types<Line>::Index;
using Size = Types<Line>::Size;
using Vector = Types<Line>::Vector;
using Domain = Types<Line>::Domain;
using BoundaryVals = Types<Line>::BoundaryVals;
using Function = Types<Line>::Function;
using PDEParams = Types<Line>::PDEParams;
using SchwarzParams = Types<Line>::SchwarzParams;
using SolverParams = Types<Line>::SolverParams;
using StatusCode = Types<Line>::StatusCode;



// Utility to calculate L2 Error Norm
Real compute_L2_error(const Vector& numerical, const std::function<Real(Real)>& exact, Real h, Size N) {
    Real error_sq = 0.0;
    // Iterate over internal nodes (assuming Dirichlet boundaries are fixed/handled)
    // Adjust loop bounds based on how 'Vector' stores nodes (e.g., including or excluding boundaries)
    for (size_t i = 0; i < numerical.size(); ++i) {
        Real x = i * h; // Assuming node 0 is at x=0. Adjust if numerical solution excludes boundaries.
        Real diff = numerical[i] - exact(x);
        error_sq += diff * diff;
    }
    return std::sqrt(error_sq * h);
}

int main() {
    std::cout << "===========================================" << std::endl;
    std::cout << "  Overlapping Schwarz Solver Test (1D)     " << std::endl;
    std::cout << "===========================================" << std::endl;

    // -----------------------------------------------------
    // 1. Problem Setup: -u''(x) = f(x) on [0,1]
    // -----------------------------------------------------
    // Exact solution: u(x) = sin(pi * x)
    // Derivatives:    u'(x) = pi * cos(pi * x)
    //                 u''(x) = -pi^2 * sin(pi * x)
    // PDE:            -u'' = pi^2 * sin(pi * x) -> f(x)
    // Parameters:     mu = 1.0, c = 0.0
    // -----------------------------------------------------
    
    Real mu = 1.0;
    Real c = 0.0;
    
    // Domain [0, 5]
    Domain omega = {0.0, 5.0}; 
    
    // Boundary Values u(0)=0, u(1)=0
    BoundaryVals dirichlet_bcs = {0.0, 0.0};

    // Forcing function f(x)
    Function forcing_function = [](Real x) {
        return M_PI * M_PI * std::sin(M_PI * x);
    };

    // Exact solution for validation
    auto exact_solution = [](Real x) {
        return std::sin(M_PI * x);
    };

    // -----------------------------------------------------
    // 2. Discretization & Solver Parameters
    // -----------------------------------------------------
    Size N_nodes = 501;
    Real h = (omega.b - omega.a) / (N_nodes - 1);
    
    // Configure PDE Parameters
    // Note: Assuming struct construction order based on typical usage.
    // You may need to adjust field names to match types.hpp definitions exactly.
    PDEParams pde_params{mu,c,forcing_function,omega,dirichlet_bcs};

    // Configure Schwarz Parameters
    SchwarzParams schwarz_params;
    schwarz_params.N = 4;           // 4 Subdomains
    schwarz_params.delta = h * 4.0;    // Overlap width (e.g., 4 grid points)

    // Configure Solver Parameters
    SolverParams solver_params;
    solver_params.max_iter = 10000;
    solver_params.eps = 1e-6;

    // -----------------------------------------------------
    // 3. Instantiate and Run Solver
    // -----------------------------------------------------
    std::cout << "Initializing DiscreteSolver..." << std::endl;
    std::cout << "\tGrid size (h):   " << h << std::endl;
    std::cout << "\tNodes:           " << N_nodes << std::endl;
    std::cout << "\tSubdomains:      " << schwarz_params.N << std::endl;
    std::cout << "\tOverlap:         " << schwarz_params.delta << std::endl;
    std::cout << "\tMax Iterations:  " << solver_params.max_iter << std::endl;
    std::cout << "\tTolerance:       " << solver_params.eps << std::endl;
    std::cout << "\tThreads:         " << omp_get_max_threads() << std::endl;


    try {
        // Instantiate the solver specialized for Line (1D)
        DiscreteSolver<Line> solver(pde_params, schwarz_params, solver_params, h);

        std::cout << "Solving system..." << std::endl;
        solver.solve();

        // Retrieve results
        Vector u_num = solver.get_solution();

        // -----------------------------------------------------
        // 4. Verification
        // -----------------------------------------------------
        
        // Check if dimensions match
        if (u_num.size() != N_nodes) {
            std::cerr << "[Warning] Solution vector size (" << u_num.size() 
                      << ") differs from expected N_nodes (" << N_nodes << ")." << std::endl;
        }

        // Calculate Error
        Real l2_error = compute_L2_error(u_num, exact_solution, h, N_nodes);


        std::cout << "\n-------------------------------------------" << std::endl;
        std::cout << "Results:" << std::endl;
        std::cout << "\tStatus:           " << (solver.status.converged() ? "CONVERGED" : "NOT CONVERGED") << std::endl;
        std::cout << "\tStatus message:   " << solver.status.message << std::endl;
        std::cout << "\tIterations:       " << solver.iter << std::endl;
        std::cout << "\tResidual:         " << solver.iter_diff << std::endl;
        std::cout << "\tL2 Error Norm:    " << l2_error << std::endl;
        std::cout << "-------------------------------------------" << std::endl;


        // Assertion for automated testing
        if (l2_error < 1e-3) {
            std::cout << "[PASSED] Error is within acceptable bounds." << std::endl;
            solver.print_to_file("./outputs/solution_test_solver.vtk");
            return 0;
        } else {
            std::cerr << "[FAILED] Error is too high." << std::endl;
            return 1;
        }

    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Exception during solve: " << e.what() << std::endl;
        return 1;
    }
}