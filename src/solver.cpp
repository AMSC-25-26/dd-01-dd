#include <solver.hpp>

SubdomainSolver::SubdomainSolver(const PDEPArams &pdep, const SchwarzParams &sp, BoundaryVals &bv, const Real h, const Index i)
        : mu(pdep.mu), c(pdep.c), f(pdep.f), omega(pdep.omega),
          Nsub(sp.N), delta(sp.delta),
          boundary_values(bv),
          h(h),
          i(i)
{
    //TODO initialize FactorizedTridiag
}