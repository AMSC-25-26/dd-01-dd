#ifndef DD_01_DD_SOLVER_HPP
#define DD_01_DD_SOLVER_HPP

#include <types.hpp>

template<int dim>
class PDESolver {
protected:
    Real mu, c, eps, delta, h;
    Domain<dim> omega;
    Function f;
    int max_iter, Nsub;

protected:
    SubIndexes global_to_sub(Index k) const;
    Index sub_to_local(SubIndexes ij) const;

public:
    PDESolver(const PDEParams<dim>& pde_params,
              const SchwarzParams& schwarz_params,
              const SolverParams& solver_params);

};

#endif //DD_01_DD_SOLVER_HPP
