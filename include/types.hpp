#ifndef DD_01_DD_TYPES_HPP
#define DD_01_DD_TYPES_HPP

#include <vector>
#include <functional>
#include <cstddef>

template 
class Types {
    using Real = double;
    using Index = int;
    using Size = size_t;
    using Vector = std::vector<Real>;
    using Function = std::function<Real(Real)>;
    
    template<int dim>
    struct Domain {
        Real a, b;
    };

    template<int dim>
    struct BoundaryVals {
        Real u_a, u_b;
    };

    struct SubIndexes {
        Index i, j;
    };

    template<int dim>
    struct PDEParams {
        Real mu, c;
        Function f;
        Domain<dim> omega;
    };

    struct SchwarzParams {
        int N;
        Real delta;
    };

    struct SolverParams {
        Real eps;
        int maxiter;
    };
}