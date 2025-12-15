#ifndef DD_01_DD_TYPES_HPP
#define DD_01_DD_TYPES_HPP

#include <vector>
#include <functional>
#include <string>

#define floating_point_error_tolerance 1e-15

enum Dimension {
    Line = 1,
    Plane = 2
};

template<Dimension dim> class Types;

template<> class Types<Line> {
    public:
    using Real = double;
    using Index = int;
    using Size = int;
    using Vector = std::vector<Real>;       //TODO decide if std::vector or Eigen::VectorXd
    using Function = std::function<Real(Real)>;

    enum StatusCode {
        Ok,
        MaxIterReached,
        SolveNotAttempted
    };

    class Status {
    public:
        std::string message;
        StatusCode code;

        bool converged() const { return code == StatusCode::Ok; }

    };

    struct Domain {
        Real a, b;
    };

    struct BoundaryVals{
        Real u_a, u_b; 
    };

    struct Boundary {
        Real a, b;
    };

    struct SubIndexes {
        Index i, j;
    };

    struct PDEParams {
        Real mu, c;
        Function f;
        Domain omega;
        BoundaryVals dirichlet;
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

template<> class Types<Plane> {
    /** TODO: part two of project **/
};

#endif //DD_01_DD_TYPES_HPP