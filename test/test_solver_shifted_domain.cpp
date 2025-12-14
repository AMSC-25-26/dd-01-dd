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

static Real compute_L2_error(const Vector& numerical, const std::function<Real(Real)> &exact, Real h, Real a) {
    Real error_sq = 0.0;
    for (size_t i = 0; i < numerical.size(); ++i) {
        Real x = a + static_cast<Real>(i) * h;
        Real diff = numerical[i] - exact(x);
        error_sq += diff * diff;
    }
    return std::sqrt(error_sq * h);
}

int main() {
    std::cout << "===========================================\n";
    std::cout << "  Schwarz Solver Test: shifted domain\n";
    std::cout << "===========================================\n";

    const Real pi = 3.141592653589793;
    const Domain omega = {0.0, 2.0};
    const BoundaryVals dirichlet_bcs = {0.0, 0.0};

    Size N_nodes = 201;
    Real h = (omega.b - omega.a) / (N_nodes - 1);

    // Exact: u(x) = sin(pi x / 2) satisfies u(0)=u(2)=0
    auto exact_solution = [pi](Real x) {
        return std::sin(0.5 * pi * x);
    };

    PDEParams pde_params{};
    pde_params.mu = 1.0;
    pde_params.c = 0.0;
    pde_params.omega = omega;
    pde_params.dirichlet = dirichlet_bcs;
    pde_params.f = [pi](Real x) {
        return (pi * pi / 4.0) * std::sin(0.5 * pi * x);
    };

    SchwarzParams schwarz_params{};
    schwarz_params.N = 4;
    schwarz_params.delta = h * 6.0;

    SolverParams solver_params{};
    solver_params.max_iter = 20000;
    solver_params.eps = 1e-8;

    DiscreteSolver<Line> solver(pde_params, schwarz_params, &solver_params, h);

    std::cout << "Solving system (shifted domain) ..." << std::endl;
    solver.solve();

    Vector u_num = solver.get_solution();
    if (u_num.size() != N_nodes) {
        std::cerr << "[Warning] Solution vector size (" << u_num.size()
                  << ") differs from expected N_nodes (" << N_nodes << ")." << std::endl;
    }

    Real l2_error = compute_L2_error(u_num, exact_solution, h, omega.a);

    std::cout << "\nResults:" << std::endl;
    std::cout << "  Status: " << (solver.status.converged() ? "CONVERGED" : "NOT CONVERGED") << std::endl;
    std::cout << "  Status message: " << solver.status.message << std::endl;
    std::cout << "  Iterations:       " << solver.status.iter << std::endl;
    std::cout << "  L2 Error Norm:    " << l2_error << std::endl;

    if (l2_error < 1e-3) {
        std::cout << "[PASSED] Shifted-domain test within tolerance." << std::endl;
        return 0;
    }
    std::cerr << "[FAILED] Error is too high." << std::endl;
    return 1;
}
