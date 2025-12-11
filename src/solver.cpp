#include <solver.hpp>

PDESolver<1>::PDESolver(const PDEParams &pde_params, const SchwarzParams &schwarz_params,
                        const SolverParams &solver_params, Real h) :
            mu(pde_params.mu),c(pde_params.c),eps(solver_params.eps), delta(schwarz_params.delta),
            h(h),omega(pde_params.omega),
            f(pde_params.f),max_iter(solver_params.max_iter),
            Nsub(schwarz_params.N)
{}

/** TODO change definition of PDESolver and instantiate normally */
SubdomainSolver<1>::SubdomainSolver(const PDEParams &pdep, const SchwarzParams &sp, BoundaryVals *bv, const Real h,
                                    const Index i) : PDESolver<1>(pdep,sp, {0.0,0},h), i(i), boundary_values(bv) {
    /** TODO actually create constructor */
    ftd = nullptr;
}

DiscreteSolver<1>::DiscreteSolver(
    const PDEParams &pdep, const SchwarzParams &sp, SolverParams *solver_params, const Real h
) : PDESolver<1>(pdep, sp, *solver_params, h) {
    subdomain_solvers.reserve(Nsub);
    
    Real subdomain_area_nonoverlapping = (omega.b - omega.a) / Nsub;

    // create a vector of SubdomainSolvers
    for (auto i = 0; i < Nsub; ++i) {
        Real a_i = ((subdomain_area_nonoverlapping * i) + omega.a) - sp.delta/2;
        Real b_i = ((subdomain_area_nonoverlapping * (i+1)) + omega.a) + sp.delta/2;
        BoundaryVals bv = {a_i, b_i};
        subdomain_solvers.emplace_back(
            SubdomainSolver<1>(pdep, sp, &bv, h, i)
        );
    }
}

Vector DiscreteSolver<1>::solve() const {
    // Variables to store u^(k) for all subdomains
    std::vector<Vector> u_current(Nsub);
    std::vector<Vector> u_prev(Nsub);

    // Initialize u^(0) = u_a + (u_b - u_a) * (x - a) / (b - a)
    
    Real subdomain_area = omega.b - omega.a;

    

    for (int i = 0; i < Nsub; i++) {
        Real x_start = omega.a + i * (subdomain_area / Nsub) - (
            i==0 ? 0.0 : (delta * 0.5)
        );

        Size n_nodes = nodes_in_subdomain(i);

        // fill with linear guess
        for(int j=0; j<n_nodes; j++) {
            Real x_j = x_start + j * h;
            Real val_interpolated = boundary_values->u_a + (boundary_values->u_b - boundary_values->u_a) * (x_j - omega.a) / (omega.b - omega.a);
            u_prev[i].push_back(
                boundary_values->u_a + (boundary_values->u_b - boundary_values->u_a) * (x_j - omega.a) / (omega.b - omega.a)
            );
        }

    }

    for (int i = 0; i<max_iter; ++i) {
    }
}