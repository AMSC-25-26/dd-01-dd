#include <solver.hpp>

#include "types.hpp"

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
