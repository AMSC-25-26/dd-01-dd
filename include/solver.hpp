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
#include <tridiag_utils.hpp>
#include <filesystem>

template<Dimension dim> class PDESolver : public Types<dim>{};

/**
 * @class PDESolver
 * @brief An "abstract" class, not meant to be initialized
 * This class is not meant to be initializable and only hosts
 * functions and parameters that are shared by SubdomainSolver and DiscreteSolver
 * @see SubdomainSolver
 * @see DiscreteSolver
 */
template<> class PDESolver<Line> : public Types<Line> {

public:
    /**
     * @brief a real number strictly greater than 0.0
     * The diffusion parameter of a Elliptic PDE it must be \f$\mu > 0\f$
     */
    Real mu;
    /**
     * @brief a real number greater or equal 0.0
     * The reaction parameter of a Elliptic PDE it must be \f$c >= 0\f$
     */
    Real c;
    /**
     * @brief a real number strictly greater than 0.0
     * The overlap of the Schwarz subdomains, such that, in 1D
     * \f$b_i-a_{i+1}=\delta \foralli : 0<i<N_{subdomains}+1\f$
     */
    Real delta;
    /**
     * @brief a real number strictly greater than 0.0
     * The distance between consecutive nodes in the grid, must be strictly positive
     * and must divide the interval \f$\lefta,b\right\f$ evenly.
     */
    Real h;
    /**
     * @brief a domain on which the PDE is defined
     * In 1D is a set of 2 real numbers (a and b)
     * representing interval \f$\lefta,b\right\f$
     * @see Types
     */
    Domain omega;
    /**
     * @brief the boundary conditions of the PDE
     * In 1D is a set of 2 real numbers representing the
     * values of the solution at the boundary of omega
     * @see Types
     */
    BoundaryVals dirichlet;
    /**
     * @brief a function on the domain omega
     * A callable object. It is imposed that \f$\calLu(x) = f(x) \forall x \in \Omega\f$.
     * It is necessary for f to be defined for every x in omega
     */
    Function f;
    /**
     * @brief the number of subdomains
     * A strictly positive integer representing the number of subdomains
     * that omega will be divided into.
     */
    Size Nsub;


    PDESolver(const PDEParams& pde_params,
              const SchwarzParams& schwarz_params,
              const Real h);
    ~PDESolver() = default;


    protected:
    Size Nnodes;                    //*< @brief the number of total nodes
    Real subdomain_area;            //*< @brief area of each nonoverlapping Schwarz subdomain
    Real domain_area;               //*< @brief the total area of the domain

    /**
     * Given the index of a subdomain (starting from 0 to Nsub-1)
     * computes the position of \f$lefta_i,b_i\right\f$
     * @param i the index of the subdomain
     * @return the boundary of the selected subdomain
     */
    Boundary get_subdomain_nonoverlapping_boundary(Index i) const;

    /**
     * Given the index of a subdomain (starting from 0 to Nsub-1)
     * computes the position of \f$lefta_i-\frac{\delta}{2},b_i+\frac{\delta}{2}\right\f$
     * @param i the index of the subdomain
     * @return the boundary of the selected subdomain
     */
    Boundary get_subdomain_overlapping_boundary(Index i) const;

    /**
     * Given a boundary returns the global index of the leftmost node inside it
     * @param boundary the boundary to examine
     * @return the index of the leftmost node
     */
    Index get_leftmost_node(Boundary boundary) const;

    /**
     * Given a boundary returns the global index of the rightmost node inside it
     * @param boundary the boundary to examine
     * @return the index of the rightmost node
     */
    Index get_rightmost_node(Boundary boundary) const;
    /**
     * Given a boundary returns the number of nodes inside it
     * @param boundary the boundary to examine
     * @return the nuber of nodes in the boundary
     */
    int get_number_of_contained_nodes(Boundary boundary) const;

    /**
     * Given a set if indices i and j, compute the global index of the
     * j-th node inside the i-th subdomain.
     * @param ij the pair of indices
     * @return the global index
     */
    Index sub_to_local(SubIndexes ij) const noexcept;

    /**
     * Given a subdomain compute the number of nodes inside it
     * @param i the index of the subdomain to be examined
     * @return the number of nodes in the subdomain
     */
    int nodes_in_subdomain(Index i) const noexcept;
};


template<Dimension dim> class SubdomainSolver : public PDESolver<dim>{};

/**
 * @class SubdomainSolver
 * @brief Solve \f$A_iu_i = b_i\f$.
 * This class solves creates, updates and compute the solution to the linear system
 * associated to each Schwarz subdomain of the PDE.
 * The solution is computed by storing the Thomas factorization of the tridiagonal matrix.
 * @see FactorizedTridiag
 */
template<> class SubdomainSolver<Line> : public PDESolver<Line> {


    public:
        Index i;                            //*< @brief This index matches the thread id and the subdomain
        /**
         * @brief Contains current values at the boundary, will be updated at each iteration
         * @see SubdomainSolver::update_boundary
         */
        BoundaryVals boundary_values;
        Vector b;                           //*< @brief the right end side of the algebaric system
        FactorizedTridiag ftd;              //*< @see FactorizedTridiag

        SubdomainSolver(const PDEParams &pdep, const SchwarzParams &sp, BoundaryVals bv, const Real h, const Index i);
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

    protected:
    Size N_overlap;                 //*< @brief the number of nodes in the subdomain including overlap
    Size N_nonoverlap;              //*< @brief the number of nodes in the subdomain without overlap
};


template<Dimension dim> class DiscreteSolver : public PDESolver<dim>{};

/**
 * @class DiscreteSolver
 * @brief Solve the Overlapping Schwarz discrete PDE problem.
 * This class manages the overall solution of the Overlapping Schwarz PDE problem
 * by coordinating the subdomain solvers in parallel and iterating until convergence.
 */
template<> class DiscreteSolver<Line> : public PDESolver<Line> {

    public:
        /**
         * @brief A status representing the state of the solver
         * @see Types
         */
        Status status;
        Real iter_diff;             //*< @brief the residual computed as \f$\max_i(||u_i^{k+1}-u_i^k||)\f$
        Index iter;                 //*< @brief the iteration reached so far

        DiscreteSolver(const PDEParams &pdep, const SchwarzParams &sp, const SolverParams &solver_params, const Real h);
        ~DiscreteSolver() = default;

        /**
         * Solves $Au = b$, where $A$ is the global stiffness matrix and $b$ the global load vector,
         * by iteratively solving the subdomain problems and updating the boundary conditions.
         * @return The global solution vector
         */
        void solve();

        /**
         * @brief getter for the solution vector
         * @return the solution at the current iteration
         */
        Vector get_solution() const;

        /**
         * @brief prints the solution to a vtk file
         * Checks that the file is in the right format and prints the current solution to
         * a .vtk file readable by paraview
         * @see https://www.paraview.org/
         */
        void print_to_file(std::filesystem::path output);
        /**
         * @brief prints the solution to standard output
         * @see https://www.paraview.org/
         */
        void print() const;

    protected:
        int max_iter;                                               //*< @brief the maximum nuber of iteration before stopping
        std::vector<SubdomainSolver<Line>> subdomain_solvers;       //*< @see SubdomainSolver
        Vector u_k;                                                 //*< @brief the currently computed solution
        Vector u_next;                                              //*< @brief a vector to host the next iteration while it's being computed
        Real eps;                                                   //*< @brief the tolerance s.t. \f$||\delta^k|| < \epsilon\f$


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
        // Vector gather_partial_solutions() const;

        /**
         * Computes the current boundary conditions for subdomain i
         * @param i The subdomain index
         * @return Computed boundary conditions
         */
        BoundaryVals current_boundary_cond(Index i) const;
};

#endif //DD_01_DD_SOLVER_HPP