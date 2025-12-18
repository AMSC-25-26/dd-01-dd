#include <types.hpp>
#include <solver.hpp>

#include <cmath>
#include <iostream>

using Real = Types<Line>::Real;
using Size = Types<Line>::Size;
using Vector = Types<Line>::Vector;
using Domain = Types<Line>::Domain;
using BoundaryVals = Types<Line>::BoundaryVals;
using PDEParams = Types<Line>::PDEParams;
using SchwarzParams = Types<Line>::SchwarzParams;
using SolverParams = Types<Line>::SolverParams;

static Real compute_L2_error(const Vector& numerical, const std::function<Real(Real)> &exact, Real h) {
    Real error_sq = 0.0;
    for (size_t i = 0; i < numerical.size(); ++i) {
        Real x = static_cast<Real>(i) * h;
        Real diff = numerical[i] - exact(x);
        error_sq += diff * diff;
    }
    return std::sqrt(error_sq * h);
}

int main() {
    std::cout << "===========================================\n";
    std::cout << "  Schwarz Solver Test: cubic-like profile\n";
    std::cout << "===========================================\n";

    const Domain omega = {0.0, 1.0};
    const BoundaryVals dirichlet_bcs = {0.0, 0.0};

    Size N_nodes = 101;
    Real h = (omega.b - omega.a) / (N_nodes - 1);

    // Exact: u(x) = x^2(1 - x) -> u'' = 2(1 - 3x), so -u'' = -2 + 6x
    auto exact_solution = [](Real x) {
        return x * x * (1.0 - x);
    };

    PDEParams pde_params{};
    pde_params.mu = 1.0;
    pde_params.c = 0.0;
    pde_params.omega = omega;
    pde_params.dirichlet = dirichlet_bcs;
    pde_params.f = [](Real x) {
        return -2.0 + 6.0 * x;
    };

    SchwarzParams schwarz_params{};
    schwarz_params.N = 4;
    schwarz_params.delta = h * 4.0;

    SolverParams solver_params{};
    solver_params.max_iter = 20000;
    solver_params.eps = 1e-8;

    DiscreteSolver<Line> solver(pde_params, schwarz_params, solver_params, h);

    std::cout << "Solving system (polynomial profile) ..." << std::endl;
    solver.solve();

    Vector u_num = solver.get_solution();
    if (u_num.size() != N_nodes) {
        std::cerr << "[Warning] Solution vector size (" << u_num.size()
                  << ") differs from expected N_nodes (" << N_nodes << ")." << std::endl;
    }

    Real l2_error = compute_L2_error(u_num, exact_solution, h);

    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "Results:" << std::endl;
    std::cout << "\tStatus:           " << (solver.status.converged() ? "CONVERGED" : "NOT CONVERGED") << std::endl;
    std::cout << "\tStatus message:   " << solver.status.message << std::endl;
    std::cout << "\tIterations:       " << solver.iter << std::endl;
    std::cout << "\tResidual:         " << solver.iter_diff << std::endl;
    std::cout << "\tL2 Error Norm:    " << l2_error << std::endl;
    std::cout << "-------------------------------------------" << std::endl;

    if (l2_error < 1e-3) {
        std::cout << "[PASSED] Polynomial test within tolerance." << std::endl;
        return 0;
    }
    std::cerr << "[FAILED] Error is too high." << std::endl;
    return 1;
}
