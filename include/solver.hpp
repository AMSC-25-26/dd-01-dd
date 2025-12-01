/******************************************************************************
 * @file solver.hpp
 * @brief This header contains classes that solve the Overlapping Schwarz
 *        PDE problem
 * @see PDESolver
 * @see DiscreteSolver
 * @see SubdomainSolver
 *******************************************************************************/


#ifndef DD_01_DD_SOLVER_HPP
#define DD_01_DD_SOLVER_HPP

#include <types.hpp>
#include <TridiagUtils.hpp>

template<size_t dim> class PDESolver : protected Types<dim>{};

template<> class PDESolver<1> : protected Types<1> {
    public:
    Real mu, c, eps, delta, h;
    Domain omega;
    Function f;
    /** TODO move these to discrete solver along with SolverParams */
    int max_iter, Nsub;

    protected:
    SubIndexes global_to_sub(Index k) const;
    Index sub_to_local(SubIndexes ij) const;

    /** TODO Maybe remove constructor */
    public:
    PDESolver(const PDEParams& pde_params,
              const SchwarzParams& schwarz_params,
              const SolverParams& solver_params,
              Real h);
    ~PDESolver() = default;

};


template<size_t dim> class SubdomainSolver : protected PDESolver<dim>{};

/**
 * @class SubdomainSolver
 * @brief Solve \f$A_iu_i = b_i\f$.
 * This class solves creates, updates and compute the solution to the linear system
 * associated to each Schwarz subdomain of the PDE.
 * The solution is computed by storing the Thomas factorization of the tridiagonal matrix.
 * @see FactorizedTridiag
 */
template<> class SubdomainSolver<1> : protected PDESolver<1> {


    public:
        Index i;                            //*< @brief This index matches the thread id and the subdomain
        /**
         * @brief Contains current values at the boundary, will be updated at each iteration
         * @see SubdomainSolver::update_boundary
         */
        BoundaryVals *boundary_values{};
        FactorizedTridiag *ftd;             //*< @see FactorizedTridiag

        SubdomainSolver(const PDEParams &pdep, const SchwarzParams &sp, BoundaryVals *bv, const Real h, const Index i);
        ~SubdomainSolver() = default;

        /**
         * Solve for \f$A_iu_i = b_i\f$, where \f$b_i\f$ is a vector of length \f$1+n_i+1\f$.
         * Where \f$b_i = \left u_{i-1,n_i-1-l} f_j u_{i+1,l} \right^T\f$,
         * given that \f$f_j = f(x_{\sigma_i(j)}) \forallj<=n_i\f$
         * and \f$l\f$ is the overlap index
         *
         * @return The solution \f$u_i\f$ on the associated subdomain.
         */
        Vector solve() const;

        /**
         * @brief Actually compute the tridiagonal factorization.
         * @see https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
         */
        void factorize();

        /**
         * @brief Updates the boundary values. This needs to be called at each iteration.
         * @param bv The new boundary values
         */
        void update_boundary(BoundaryVals bv);

    private:
        /**
         * @brief Computes \f$(A_i)_{j,k}\f$
         * Computes the coefficient at position (j,k) of the local tridiagonal matrix
         * @param j Row index
         * @param k Collumn index
         * @return The coefficient at position j,k of the local stiffness matrix
         */
        Real stiff_mat(Index j, Index k) const;
};

#endif //DD_01_DD_SOLVER_HPP