#ifndef DD_01_DD_TYPES_HPP
#define DD_01_DD_TYPES_HPP

#include <vector>
#include <functional>
#include <type_traits>

template<size_t dim> class Types;

template<> class Types<1> {
    public:
    using Real = double;
    using Index = int;
    using Size = size_t;
    using Vector = std::vector<Real>;       //TODO decide if std::vector or Eigen::VectorXd
    using Function = std::function<Real(Real)>;

    struct Domain {
        Real a, b;
    } ;

    struct BoundaryVals{
        Real u_a, u_b; 
    };

    struct SubIndexes {
        Index i, j;
    };

    struct PDEParams {
        Real mu, c;
        Function f;
        Domain omega;
    };

    struct SchwarzParams {
        int N;
        Real delta;
    };

    struct SolverParams {
        Real eps;
        int max_iter;
    };
};

template<> class Types<2> {
    /** TODO: part two of project **/
};

#endif //DD_01_DD_TYPES_HPP