#include <types.hpp>
#include <solver.hpp>

#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include <cmath>
#include <functional>

using Real = Types<Line>::Real;
using Size = Types<Line>::Size;
using Domain = Types<Line>::Domain;
using BoundaryVals = Types<Line>::BoundaryVals;
using Function = Types<Line>::Function;
using PDEParams = Types<Line>::PDEParams;
using SchwarzParams = Types<Line>::SchwarzParams;
using SolverParams = Types<Line>::SolverParams;

struct ProblemConfig {
    std::string name;
    Size nodes;
    Real eps;
};

void run_benchmark(const ProblemConfig& config, const std::vector<int>& subdomain_counts) {
    Real mu = 1.0;
    Real c = 0.0;
    Domain omega = {0.0, 10.0}; 
    BoundaryVals dirichlet_bcs = {0.0, 0.0};
    
    Function forcing_function = [](Real x) -> Real {
        return M_PI * M_PI * std::sin(M_PI * x);
    };

    Real h = (omega.b - omega.a) / (config.nodes - 1);
    PDEParams pde_params{mu, c, forcing_function, omega, dirichlet_bcs};
    
    SolverParams solver_params;
    solver_params.max_iter = 1000000;
    solver_params.eps = config.eps;

    std::cout << "\n>>> BENCHMARK: " << config.name << " (" << config.nodes << " nodes)" << std::endl;
    std::cout << std::left << std::setw(12) << "Subdomains" 
              << std::setw(10) << "Threads" 
              << std::setw(12) << "Status"
              << std::setw(10) << "Iters" 
              << std::setw(15) << "Time (ms)" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    for (int N : subdomain_counts) {
        std::vector<int> test_threads = (N == 1) ? std::vector<int>{1} : std::vector<int>{1, N};

        for (int n_threads : test_threads) {
            omp_set_num_threads(n_threads);
            SchwarzParams schwarz_params;
            schwarz_params.N = N;
            schwarz_params.delta = h * 4.0; 

            try {
                DiscreteSolver<Line> solver(pde_params, schwarz_params, solver_params, h);

                auto start = std::chrono::high_resolution_clock::now();
                solver.solve();
                auto end = std::chrono::high_resolution_clock::now();

                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                std::string status_str = solver.status.converged() ? "CONVERGED" : "FAILED";

                std::cout << std::left << std::setw(12) << N 
                          << std::setw(10) << n_threads 
                          << std::setw(12) << status_str
                          << std::setw(10) << solver.iter 
                          << std::setw(15) << duration << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << std::endl;
            }
        }
    }
}

int main() {
    std::cout << "==============================================================================" << std::endl;
    std::cout << "  Schwarz Solver: Comparative Benchmarking (Base vs Heavy) " << std::endl;
    std::cout << "==============================================================================" << std::endl;

    std::vector<int> subdomains = {1, 2, 4, 8, 16};

    // Base Problem: Similar to your previous run
    ProblemConfig base_problem = {"Base Problem", 5001, 1e-4};
    
    // Heavy Problem: 10x more nodes (50k nodes)
    // This significantly increases the local work per subdomain and memory overhead
    ProblemConfig heavy_problem = {"Heavy Problem", 50001, 1e-4};

    run_benchmark(base_problem, subdomains);
    run_benchmark(heavy_problem, subdomains);

    std::cout << "\n==============================================================================" << std::endl;
    return 0;
}