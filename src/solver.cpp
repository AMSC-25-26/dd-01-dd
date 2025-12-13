#include <solver.hpp>
#include <types.hpp>

Types<Line>::Boundary PDESolver<Line>::get_subdomain_nonoverlapping_boundary(Index i) const {
    return {((subdomain_area * i) + omega.a),((subdomain_area * (i+1)) + omega.a)};
}

Types<Line>::Boundary PDESolver<Line>::get_subdomain_overlapping_boundary(Index i) const {
    Boundary boundary = get_subdomain_nonoverlapping_boundary(i);
    if (i!=0) boundary.a -= delta*.5;
    if (i!=Nsub-1) boundary.b += delta*.5;
    return boundary;
}

Types<Line>::Index PDESolver<Line>::get_leftmost_node(Boundary boundary) const {
    return static_cast<Index>((boundary.a-omega.a)/h)+1;
}

Types<Line>::Index PDESolver<Line>::get_rightmost_node(Boundary boundary) const {
    return static_cast<Index>((boundary.b-omega.a)/h);
}

Types<Line>::Index PDESolver<Line>::get_number_of_contained_nodes(Boundary boundary) const {
    return get_rightmost_node(boundary)-get_leftmost_node(boundary);
}

PDESolver<Line>::PDESolver(const PDEParams &pde_params, const SchwarzParams &schwarz_params, Real h) : mu(pde_params.mu),
                                                                                                       c(pde_params.c), delta(schwarz_params.delta),
                                                                                                       h(h), omega(pde_params.omega), dirichlet(pde_params.dirichlet),
                                                                                                       f(pde_params.f),
                                                                                                       Nsub(schwarz_params.N) {
    domain_area = omega.b - omega.a;
    subdomain_area = domain_area/Nsub;
    Nnodes = static_cast<int>((omega.b - omega.a) / h) + 1;
    /* TODO check conditions:
     * b-a is perfectly divisible by h
     * need b-a != 0
     * maybe more...
     */
}

/** TODO change definition of PDESolver and instantiate normally */
SubdomainSolver<Line>::SubdomainSolver(const PDEParams &pdep, const SchwarzParams &sp, BoundaryVals bv, const Real h,
                                    const Index i) : PDESolver<Line>(pdep, sp, h), i(i), boundary_values(bv), ftd(0) {
    N_overlap = get_number_of_contained_nodes(
        get_subdomain_overlapping_boundary(i)
        );
    N_nonoverlap = get_number_of_contained_nodes(
        get_subdomain_nonoverlapping_boundary(i)
        );
    b.resize(N_overlap);
    ftd = FactorizedTridiag(N_overlap);

    ftd(0,0) = 1;
    ftd(N_overlap-1,N_overlap-1) = 1;
    for (Index j = 1; j < (Index)N_overlap-1; ++j) {
        ftd(j,j-1) = -mu/(h*h);
        ftd(j,j) = (2*mu/(h*h))+c;
        ftd(j,j+1) = -mu/(h*h);
    }

    b[0] = bv.u_a;
    b[N_overlap-1] = bv.u_b;
    for (Index j = 1; j < (Index)N_overlap-1; ++j)
        b[j] = f(this->sub_to_local({i,j}) * h + this->omega.a);
}

Types<Line>::Vector SubdomainSolver<Line>::solve() const {
    return ftd.solve(b);
}

void SubdomainSolver<Line>::factorize() {
    ftd.factorize();
}

void SubdomainSolver<Line>::update_boundary(BoundaryVals bv) {
    boundary_values = bv;
    b[0] = bv.u_a;
    b[N_overlap-1] = bv.u_b;
}

DiscreteSolver<Line>::DiscreteSolver(
    const PDEParams &pdep, const SchwarzParams &sp, SolverParams *solver_params, const Real h
) : PDESolver<Line>(pdep, sp, h), max_iter(solver_params->max_iter), iter(0), iter_diff(0), eps(solver_params->eps) {

    status.code = SolveNotAttempted;
    status.message = "You have yet to call solve()";

    u_k.resize(Nnodes);
    u_next.resize(Nnodes);
    subdomain_solvers.reserve(Nsub);

    // create a vector of SubdomainSolvers
    BoundaryVals bv = {0.0,0.0};
    for (auto i = 0; i < Nsub; ++i) {
        subdomain_solvers.emplace_back(pdep, sp, bv, h, i);
    }
}

void DiscreteSolver<Line>::solve() {
    // Initialize u^(0) = u_a + (u_b - u_a) * (x - a) / (b - a)
    Real slope = (dirichlet.u_b-dirichlet.u_a)/(omega.b-omega.a);
    for (Size i = 0; i < Nnodes; ++i) {
        u_k[i] = (static_cast<Real>(i)*h - omega.a)*slope + dirichlet.u_a;
    }

    iter_diff = eps + 1;
    while (iter++ < max_iter && iter_diff > eps) {
        advance();
    }

    if (max_iter == iter) {
        status.message = "Maximum number of iterations exceeded";
        status.code = MaxIterReached;
    } else {
        status.message = "Solver reached solution without problems";
        status.code = Ok;
    }
}

Types<Line>::Vector DiscreteSolver<Line>::get_solution() const {
    return u_k;
}

void DiscreteSolver<Line>::advance() {
    u_next.clear();
    for(auto i=0; i<Nsub; ++i) {
        if (iter == 0) {
            subdomain_solvers[i].factorize();
        }
        subdomain_solvers[i].update_boundary(
            current_boundary_cond(i)
        );
        Vector u_i_k = subdomain_solvers[i].solve();
        u_next.insert(u_next.end(), u_i_k.begin(), u_i_k.end());
    }
    std::swap(u_k, u_next);
}

Types<Line>::BoundaryVals DiscreteSolver<Line>::current_boundary_cond(Index i) const {
    const Boundary boundary = get_subdomain_overlapping_boundary(i);
    if (i == 0) {
        return {dirichlet.u_a,u_k[get_rightmost_node(boundary)]};
    }
    if (i == Nsub-1) {
        return {u_k[get_leftmost_node(boundary)],dirichlet.u_b};
    }
    return {
        u_k[get_leftmost_node(boundary)],
        u_k[get_rightmost_node(boundary)]
    };
}

