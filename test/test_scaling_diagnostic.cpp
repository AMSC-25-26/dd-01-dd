#include <types.hpp>
#include <solver.hpp>

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

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

const Real PI = 3.141592653589793;

Real compute_L2_error(const Vector& numerical, const std::function<Real(Real)>& exact, Real h, Real omega_a) {
    Real error_sq = 0.0;
    for (size_t i = 0; i < numerical.size(); ++i) {
        Real x = omega_a + i * h;
        Real diff = numerical[i] - exact(x);
        error_sq += diff * diff;
    }
    return std::sqrt(error_sq * h);
}

Real compute_max_error(const Vector& numerical, const std::function<Real(Real)>& exact, Real h, Real omega_a) {
    Real max_err = 0.0;
    for (size_t i = 0; i < numerical.size(); ++i) {
        Real x = omega_a + i * h;
        Real err = std::abs(numerical[i] - exact(x));
        if (err > max_err) max_err = err;
    }
    return max_err;
}

struct TestResult {
    Size N_nodes;
    int Nsub;
    Real delta;
    Real h;
    Real L2_error;
    Real max_error;
    int iterations;
    bool converged;
    Size computed_Nnodes;
    std::string notes;
};

TestResult run_test(Size N_nodes, int Nsub, int overlap_nodes, Real domain_length = 1.0) {
    TestResult result;
    result.N_nodes = N_nodes;
    result.Nsub = Nsub;

    Real h = domain_length / (N_nodes - 1);
    result.h = h;
    result.delta = h * overlap_nodes;

    Domain omega = {0.0, domain_length};
    BoundaryVals dirichlet_bcs = {0.0, 0.0};

    Function forcing_function = [](Real x) {
        return PI * PI * std::sin(PI * x);
    };

    auto exact_solution = [](Real x) {
        return std::sin(PI * x);
    };

    PDEParams pde_params;
    pde_params.mu = 1.0;
    pde_params.c = 0.0;
    pde_params.omega = omega;
    pde_params.dirichlet = dirichlet_bcs;
    pde_params.f = forcing_function;

    SchwarzParams schwarz_params;
    schwarz_params.N = Nsub;
    schwarz_params.delta = result.delta;

    SolverParams solver_params;
    solver_params.max_iter = 10000;
    solver_params.eps = 1e-8;

    try {
        DiscreteSolver<Line> solver(pde_params, schwarz_params, &solver_params, h);
        solver.solve();

        Vector u_num = solver.get_solution();
        result.computed_Nnodes = u_num.size();
        result.L2_error = compute_L2_error(u_num, exact_solution, h, omega.a);
        result.max_error = compute_max_error(u_num, exact_solution, h, omega.a);
        result.iterations = solver.status.iter;
        result.converged = solver.status.converged();

        // Check for issues
        if (result.computed_Nnodes != N_nodes) {
            result.notes = "SIZE MISMATCH: expected " + std::to_string(N_nodes) +
                          " got " + std::to_string(result.computed_Nnodes);
        }

        // Check divisibility
        Size nodes_per_subdomain = (N_nodes - 1) / Nsub;
        Size remainder = (N_nodes - 1) % Nsub;
        if (remainder != 0) {
            result.notes += " NOT_DIVISIBLE: (N-1)%" + std::to_string(Nsub) + "=" + std::to_string(remainder);
        }

    } catch (const std::exception& e) {
        result.L2_error = -1;
        result.max_error = -1;
        result.iterations = -1;
        result.converged = false;
        result.notes = std::string("EXCEPTION: ") + e.what();
    }

    return result;
}

void print_result(const TestResult& r) {
    std::cout << std::setw(8) << r.N_nodes
              << std::setw(6) << r.Nsub
              << std::setw(10) << std::scientific << std::setprecision(2) << r.h
              << std::setw(10) << r.delta
              << std::setw(12) << r.L2_error
              << std::setw(12) << r.max_error
              << std::setw(8) << std::defaultfloat << r.iterations
              << std::setw(10) << (r.converged ? "YES" : "NO")
              << std::setw(8) << r.computed_Nnodes
              << "  " << r.notes
              << std::endl;
}

int main() {
    std::cout << "==============================================================================" << std::endl;
    std::cout << "  Scaling Diagnostic Test for Overlapping Schwarz Solver" << std::endl;
    std::cout << "  Domain: [0, 1], PDE: -u'' = pi^2*sin(pi*x), exact: u(x) = sin(pi*x)" << std::endl;
    std::cout << "==============================================================================" << std::endl;

    std::cout << "\n--- Header ---" << std::endl;
    std::cout << std::setw(8) << "N_nodes"
              << std::setw(6) << "Nsub"
              << std::setw(10) << "h"
              << std::setw(10) << "delta"
              << std::setw(12) << "L2_err"
              << std::setw(12) << "max_err"
              << std::setw(8) << "iter"
              << std::setw(10) << "conv"
              << std::setw(8) << "size"
              << "  notes"
              << std::endl;
    std::cout << std::string(100, '-') << std::endl;

    // Test 1: Baseline - 1 subdomain (no decomposition)
    std::cout << "\n[Test 1: Single subdomain - baseline]" << std::endl;
    for (Size N : {11, 21, 51, 101, 201, 501, 1001}) {
        print_result(run_test(N, 1, 0));
    }

    // Test 2: 2 subdomains - check divisibility
    std::cout << "\n[Test 2: 2 subdomains - divisibility check]" << std::endl;
    for (Size N : {11, 21, 31, 41, 51, 101, 201, 501}) {
        print_result(run_test(N, 2, 4));
    }

    // Test 3: 4 subdomains
    std::cout << "\n[Test 3: 4 subdomains]" << std::endl;
    for (Size N : {21, 41, 81, 101, 201, 401, 501, 1001}) {
        print_result(run_test(N, 4, 4));
    }

    // Test 4: Many subdomains
    std::cout << "\n[Test 4: 8 subdomains]" << std::endl;
    for (Size N : {41, 81, 161, 321, 801, 1001}) {
        print_result(run_test(N, 8, 4));
    }

    // Test 5: Specific problematic cases - N_nodes that are NOT (k*Nsub + 1)
    std::cout << "\n[Test 5: Problematic N_nodes values with 4 subdomains]" << std::endl;
    std::cout << "  (N-1) should be divisible by Nsub for clean decomposition" << std::endl;
    // Good: N = k*4 + 1 = 5, 9, 13, 17, 21, 25...
    // Bad: N = 6, 7, 8, 10, 11, 12...
    for (Size N : {5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 21, 25, 100, 101, 102, 103, 104, 105}) {
        print_result(run_test(N, 4, 2));
    }

    // Test 6: Overlap size impact
    std::cout << "\n[Test 6: Overlap size impact (N=101, Nsub=4)]" << std::endl;
    for (int overlap : {1, 2, 4, 8, 10, 15, 20}) {
        print_result(run_test(101, 4, overlap));
    }

    // Test 7: Large scale
    std::cout << "\n[Test 7: Large scale tests]" << std::endl;
    print_result(run_test(1001, 4, 4));
    print_result(run_test(1001, 8, 4));
    print_result(run_test(1001, 10, 4));
    print_result(run_test(2001, 4, 4));
    print_result(run_test(2001, 8, 4));
    print_result(run_test(5001, 4, 4));
    print_result(run_test(5001, 10, 4));
    return 0;
}


