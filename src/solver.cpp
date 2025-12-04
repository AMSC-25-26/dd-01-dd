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
    subdomain_solvers.reserve(sp.N);
    
    Real subdomain_dim_nonoverlapping = (omega.b - omega.a) / sp.N;

    // create a vector of SubdomainSolvers
    for (auto i = 0; i < sp.N; ++i) {
        Real a_i = ((subdomain_dim_nonoverlapping * i) + omega.a) - sp.delta/2;
        Real b_i = ((subdomain_dim_nonoverlapping * (i+1)) + omega.a) + sp.delta/2;
        BoundaryVals bv = {a_i, b_i};
        subdomain_solvers.emplace_back(
            SubdomainSolver<1>(pdep, sp, &bv, h, i)
        );
    }
}

Vector DiscreteSolver<1>::solve() const {
    for (int i = 0; i<max_iter; ++i) {
        for (auto &subsolver : subdomain_solvers) {
            subsolver.factorize();
            subsolver.solve();
            subsolver.update_boundary();
        }
    }
}