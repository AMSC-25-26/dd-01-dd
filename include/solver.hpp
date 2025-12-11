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

template<Dimension dim> class PDESolver : protected Types<dim>{};

template<> class PDESolver<Line> : protected Types<Line> {
    public:
    Real mu, c, eps, delta, h;
    Domain omega;
    BoundaryVals dirichlet;
    Function f;
    /** TODO move these to discrete solver along with SolverParams */
    int max_iter;
    Size Nsub, Nnodes;

    protected:
    SubIndexes global_to_sub(Index k) const noexcept;
    Index sub_to_local(SubIndexes ij) const noexcept;
    Size nodes_in_subdomain(Index i) const noexcept;

    /** TODO Maybe remove constructor */
    public:
    PDESolver(const PDEParams& pde_params,
              const SchwarzParams& schwarz_params,
              const SolverParams& solver_params,
              Real h);
    ~PDESolver() = default;

};


template<Dimension dim> class SubdomainSolver : protected PDESolver<dim>{};

/**
 * @class SubdomainSolver
 * @brief Solve \f$A_iu_i = b_i\f$.
 * This class solves creates, updates and compute the solution to the linear system
 * associated to each Schwarz subdomain of the PDE.
 * The solution is computed by storing the Thomas factorization of the tridiagonal matrix.
 * @see FactorizedTridiag
 */
template<> class SubdomainSolver<Line> : protected PDESolver<Line> {


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


template<Dimension dim> class DiscreteSolver : protected PDESolver<dim>{};

/**
 * @class DiscreteSolver
 * @brief Solve the Overlapping Schwarz discrete PDE problem.
 * This class manages the overall solution of the Overlapping Schwarz PDE problem
 * by coordinating the subdomain solvers in parallel and iterating until convergence.
 */
template<> class DiscreteSolver<Line> : protected PDESolver<Line> {

    public:
        Status status;

        DiscreteSolver(const PDEParams &pdep, const SchwarzParams &sp, SolverParams *solver_params, const Real h);
        ~DiscreteSolver() = default;

        /**
         * Solves $Au = b$, where $A$ is the global stiffness matrix and $b$ the global load vector,
         * by iteratively solving the subdomain problems and updating the boundary conditions.
         * @return The global solution vector
         */
        void solve();

        Vector get_solution() const;

        /**
         * @brief TODO
         */
        void print_to_file();

    protected:
        std::vector<SubdomainSolver<1>> subdomain_solvers;
        Vector u_k, u_next;
        Index iter;
        Real iter_diff;


    private:
        /**
         * Advances the computation of one iteration, using SubdomainSolver, on each parallel process.
         * When all processes have completed their iteration, the partial solutions are gathered.
         */
        void advance();

        /**
         * Gathers the partial solutions from each subdomain solver into a single global solution vector and
         * constucts $u_{k+1}$
         * @return The global solution vector
         */
        Vector gather_partial_solutions() const;

        /**
         * Computes the current boundary conditions for subdomain i
         * @param i The subdomain index
         * @return Computed boundary conditions
         */
        BoundaryVals current_boundary_cond(Index i) const;
};

#endif //DD_01_DD_SOLVER_HPP