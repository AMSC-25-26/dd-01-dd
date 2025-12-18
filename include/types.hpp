/******************************************************************************
 * @file types.hpp
 * @brief This header contains type definitions and structures
 * @see PDESolver
 * @see DiscreteSolver
 * @see SubdomainSolver
 *******************************************************************************/

#ifndef DD_01_DD_TYPES_HPP
#define DD_01_DD_TYPES_HPP

#include <vector>
#include <functional>
#include <string>

/** 
 * @brief Tolerance for floating point comparisons
 */
#define floating_point_error_tolerance 1e-15

/** 
 * @brief The dimension of the problem
 */
enum Dimension {
    Line = 1, ///< 1D problem
    Plane = 2 ///< 2D problem
};

/**
 * @brief A template class to host type definitions and structures
 * @tparam dim The dimension of the problem (Line or Plane)
 */
template<Dimension dim> class Types;

/**
 * @class Types
 * @brief A template class to host type definitions and structures for 1D problems
 * It defines the scalar types, vectors and parameters structures used to solve
 * \f$-\mu u'' + cu = f\f$ on an interval
 */
template<> class Types<Line> {
    public:
    using Real = double;                        ///< @brief Floating point type used for real numbers
    using Index = int;                          ///< @brief Integer type used for indices   
    using Size = int;                           ///< @brief Integer type used for sizes and counts                        
    using Vector = std::vector<Real>;           ///< @brief Structure to represent vectors of real numbers      
    using Function = std::function<Real(Real)>; ///< @brief Type representing a mathematical function \f$ f: \mathbb{R} \to \mathbb{R} \f$

    /**
     * @brief Status codes for the iterative solver
     */
    enum StatusCode {
        Ok,                 ///< The solver converged successfully
        MaxIterReached,     ///< The solver reached the maximum number of iterations without converging
        SolveNotAttempted   ///< The solver has not been run yet
    };

    /**
     * @class Status
     * @brief A structure representing the status of the solver
     */

    class Status {
    public:
        std::string message; ///< @brief A describitive message about the solver status
        StatusCode code;     ///< @brief The status code

        /**
         * @brief Checks if the solver converged successfully
         * @return true if the solver converged, false otherwise
         */
        bool converged() const { return code == StatusCode::Ok; }

    };

    /**
     * @brief Structure that represents the physical domain \f$ \Omega = (a, b) \f$
     */

    struct Domain {
        Real a, b;
    };

    /**
     * @brief Structure that represents the Dirichlet boundary conditions
     * Represents the value of the function at the boundaries of the domain
     * \f$ u(a) = u_a, u(b) = u_b \f$
     */
    struct BoundaryVals{
        Real u_a, u_b; 
    };

    /**
     * @brief Structure that represents a boundary of a subdomain
     */
    struct Boundary {
        Real a, b;
    };

    /**
     * @brief Structure that represents a pair of indices
     * This usually refers to \f$ (a_i, b_i) \f$ the boundaries of a subdomain
     */
    struct SubIndexes {
        Index i; ///< @brief The index of the subdomain
        Index j; ///< @brief The local node index within the subdomain
    };

    /**
     * @brief Structure that holds the parameters for the PDE problem
     * Encapsulates the coefficients for the equation \f$ -\mu u'' + cu = f \f$
     */
    struct PDEParams {
        Real mu;                ///< @brief Diffusion coefficient \f$ \mu > 0 \f$
        Real c;                 ///< @brief Reaction coefficient \f$ c \ge 0 \f$
        Function f;             ///< @brief Forcing term function \f$ f(x) \f$
        Domain omega;           ///< @brief The global domain \f$ \Omega \f$
        BoundaryVals dirichlet; ///< @brief The global boundary conditions
    };

    /**
     * @brief Parameters for the Schwarz domain decomposition method
     */
    struct SchwarzParams {
        int N;      ///< @brief Number of subdomains to divide \f$ \Omega \f$ into
        Real delta; ///< @brief The overlap size \f$ \delta \f$ between adjacent subdomains
    };

    /**
     * @brief Parameters for the iterative solver control
     */
    struct SolverParams {
        Real eps;     ///< @brief Tolerance \f$ \epsilon \f$ for the stopping criterion
        int max_iter; ///< @brief Maximum number of iterations allowed
    };
};

/**
 * @class Types<Plane>
 * @brief Type definitions for 2D problems (Future implementation).
 */
template<> class Types<Plane> {
    /** TODO: part two of project **/
};

#endif //DD_01_DD_TYPES_HPP