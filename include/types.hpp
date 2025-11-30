#ifndef DD_01_DD_TYPES_HPP
#define DD_01_DD_TYPES_HPP

#include <vector>
#include <functional>
#include <cstddef>

template<int dim>
class Types {
    public:
    using Real = double;
    using Index = int;
    using Size = size_t;
    using Vector = std::vector<Real>;       //TODO decide if std::vector or Eigen::VectorXd
    using Function = std::function<Real(Real)>;
 
    typedef struct {
        Real a, b; 
    } Domain1D;
    using Domain = Domain1D;

    typedef struct {
        Real u_a, u_b; 
    } BoundaryVals1D;
    using BoundaryVals = BoundaryVals1D;

    typedef struct {Index i, j; } SubIndexes;

    typedef struct {
        Real mu, c;
        Function f;
        Domain omega;
    } PDEParams;

    typedef struct {
        int N;
        Real delta;
    } SchwarzParams;

    typedef struct {
        Real eps;
        int maxiter;
    } SolverParams;
}

#endif //DD_01_DD_TYPES_HPP