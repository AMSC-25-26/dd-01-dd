#include <solver.hpp>

#include "types.hpp"

PDESolver<Line>::PDESolver(const PDEParams &pde_params, const SchwarzParams &schwarz_params,
                        const SolverParams &solver_params, Real h) :
            mu(pde_params.mu),c(pde_params.c),eps(solver_params.eps), delta(schwarz_params.delta),
            h(h),omega(pde_params.omega),
            f(pde_params.f),max_iter(solver_params.max_iter),
            Nsub(schwarz_params.N) {
    Nnodes = static_cast<int>( (omega.b-omega.a) / h )+1;
    /* TODO check conditions:
     * b-a is perfectly divisible by h
     * need b-a != 0
     * maybe more...
     */

}

/** TODO change definition of PDESolver and instantiate normally */
SubdomainSolver<Line>::SubdomainSolver(const PDEParams &pdep, const SchwarzParams &sp, BoundaryVals *bv, const Real h,
                                    const Index i) : PDESolver<Line>(pdep,sp, {0.0,0},h), i(i), boundary_values(bv) {
    /** TODO actually create constructor */
    ftd = nullptr;
}

DiscreteSolver<1>::DiscreteSolver(
    const PDEParams &pdep, const SchwarzParams &sp, SolverParams *solver_params, const Real h
) : PDESolver<1>(pdep, sp, *solver_params, h), iter(0), iter_diff(0){

    // useful renames
    u_k.reserve(Nnodes);
    u_next.reserve(Nnodes);
    subdomain_solvers.reserve(Nsub);
    
    Real subdomain_area_nonoverlapping = (omega.b - omega.a) / Nsub;

    // create a vector of SubdomainSolvers
    for (auto i = 0; i < Nsub; ++i) {
        BoundaryVals bv = {0.0, 0.0};
        subdomain_solvers.emplace_back(
            SubdomainSolver<1>(pdep, sp, &bv, h, i)
        );
    }
}

void DiscreteSolver<1>::solve() {
    // Initialize u^(0) = u_a + (u_b - u_a) * (x - a) / (b - a)
    Real slope = (dirichlet.u_b-dirichlet.u_a)/(omega.b-omega.a);
    for (Size i = 0; i < Nnodes; ++i) {
        u_k[i] = (static_cast<Real>(i)*h - omega.a)*slope + dirichlet.u_a;
    }

    iter_diff = eps + 1;
    while (iter++ < max_iter && iter_diff > eps) {
        u_next = advance();
    }

    if (max_iter == iter) {
        status.message = "Maximum number of iterations exceeded";
        status.code = MaxIterReached;
    } else {
        status.message = "Solver reached solution without problems";
        status.code = Ok;
    }
}

Types<1>::Vector DiscreteSolver<1>::get_solution() const {
    return u_k;
}

void DiscreteSolver<1>::advance() {
    u_next.clear();
    for(auto i=0; i<Nsub; ++ii) {
        if (iter == 0) {
            subdomain_solvers[i].factorize();
        }
        subdomain_solvers[i].update_boundary(
            current_boundary_cond(i)
        );
        Vector u_i_k = subdomain_solvers[i].solve();
        u_next.insert(u_next.end(), u_i_k.begin(), u_i_k.end());
    }
    std::swap(u_k, u_next);
}

Types<1>::BoundaryVals current_boundary_cond(Index i) const {
    Real a_i = ((subdomain_area_nonoverlapping * i) + omega.a) - sp.delta/2;
    Real b_i = ((subdomain_area_nonoverlapping * (i+1)) + omega.a) + sp.delta/2;
    
    // the index of the leftmost and rightmost node contained in the overlapping subdomain
    Index leftmost = static_cast<Index>(a_i/h)+1;
    Index rightmost = static_cast<Index>(b_i/h);

    return {
        u_k[leftmost],
        u_k[rightmost]
    };
}

