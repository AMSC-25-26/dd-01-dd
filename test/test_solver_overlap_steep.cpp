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
    std::cout << "  Schwarz Solver Test: steep overlap case\n";
    std::cout << "===========================================\n";

    const Domain omega = {0.0, 1.0};
    const BoundaryVals dirichlet_bcs = {0.0, 0.0};

    Size N_nodes = 151;
    Real h = (omega.b - omega.a) / (N_nodes - 1);

    PDEParams pde_params{};
    pde_params.mu = 0.5;   // lower diffusion
    pde_params.c = 5.0;    // stronger reaction
    pde_params.omega = omega;
    pde_params.dirichlet = dirichlet_bcs;
    pde_params.f = [](Real x) {
        return 0.5 * (4.0 * M_PI * M_PI) * std::sin(2.0 * M_PI * x) + 5.0 * std::sin(2.0 * M_PI * x);
    }; // -0.5 u'' + 5 u with u = sin(2 pi x)

    SchwarzParams schwarz_params{};
    schwarz_params.N = 4;
    schwarz_params.delta = h * 8.0; // generous overlap for the steep mode

    SolverParams solver_params{};
    solver_params.max_iter = 30000;
    solver_params.eps = 1e-8;

    auto exact_solution = [](Real x) {
        return std::sin(2.0 * M_PI * x);
    };

    DiscreteSolver<Line> solver(pde_params, schwarz_params, solver_params, h);

    std::cout << "Solving system (steep overlap case) ..." << std::endl;
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

    if (l2_error < 3e-3) {
        std::cout << "[PASSED] Steep test within tolerance." << std::endl;
        return 0;
    }
    std::cerr << "[FAILED] Error is too high." << std::endl;
    return 1;
}
