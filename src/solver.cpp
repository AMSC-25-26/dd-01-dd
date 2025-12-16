#include <cmath>
#include <filesystem>
#include <solver.hpp>
#include <types.hpp>

#include <fstream>
#include <iomanip>
#include <stdexcept>

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
    // check to not overshoot below 0
    if (fabs(boundary.a-omega.a) < floating_point_error_tolerance) return 0;
    return static_cast<Index>((boundary.a-omega.a)/h)+1;
}

Types<Line>::Index PDESolver<Line>::get_rightmost_node(Boundary boundary) const {
    // no need to check for overshoot
    return static_cast<Index>((boundary.b-omega.a)/h);
}

int PDESolver<Line>::get_number_of_contained_nodes(Boundary boundary) const {
    // +1 because both endpoints are included in the discrete grid
    return get_rightmost_node(boundary) - get_leftmost_node(boundary) + 1;
}

PDESolver<Line>::PDESolver(const PDEParams &pde_params, const SchwarzParams &schwarz_params, const Real h) : mu(pde_params.mu),
    c(pde_params.c), delta(schwarz_params.delta),
    h(h), omega(pde_params.omega), dirichlet(pde_params.dirichlet),
    f(pde_params.f),
    Nsub(schwarz_params.N) {
    domain_area = omega.b - omega.a;
    subdomain_area = domain_area/Nsub;
    Nnodes = static_cast<int>((omega.b - omega.a) / h) + 1;

    //checking conditions on parameters
    if (mu <= 0) throw std::invalid_argument("The diffusion parameter mu must be s.t. mu > 0");
    if (c < 0) throw std::invalid_argument("The reaction parameter c must be s.t. c >= 0 ");
    if (delta <= 0) throw std::invalid_argument("The overlap must represent a real length > 0");
    if (h <= 0) throw std::invalid_argument("The distance between nodes must be strictly positive");
    if (fabs(static_cast<Size>(domain_area/h)*h - domain_area) > floating_point_error_tolerance)
        throw std::invalid_argument("The distance between nodes must divide the domain evenly");
    if (Nsub <= 0) throw std::invalid_argument("The Schwarz method require at least 1 subdomain");


}

// get global node index from (subdomain index, local node index)rr
// TODO: verify correctness
Types<Line>::Index PDESolver<Line>::sub_to_local(SubIndexes sub) const noexcept {
    return get_leftmost_node(get_subdomain_overlapping_boundary(sub.i)) + sub.j;
}

/** TODO change definition of PDESolver and instantiate normally */
SubdomainSolver<Line>::SubdomainSolver(const PDEParams &pdep, const SchwarzParams &sp, BoundaryVals bv,const Real h,
                                       const Index i) : PDESolver<Line>(pdep, sp, h), i(i), boundary_values(bv), ftd(get_number_of_contained_nodes(get_subdomain_overlapping_boundary(i))) {
    N_overlap = get_number_of_contained_nodes(get_subdomain_overlapping_boundary(i));
    N_nonoverlap = get_number_of_contained_nodes(get_subdomain_nonoverlapping_boundary(i));
    b.resize(N_overlap);

    ftd(0,0) = 1;
    ftd(N_overlap-1,N_overlap-1) = 1;
    for (auto j = 1; j < N_overlap-1; ++j) {
        ftd(j,j-1) = -mu/(h*h);
        ftd(j,j) = (2*mu/(h*h))+c;
        ftd(j,j+1) = -mu/(h*h);
    }

    b[0] = bv.u_a;
    b[N_overlap-1] = bv.u_b;
    for (auto j = 1; j < N_overlap-1; ++j) {
        b[j] = f((sub_to_local({i,j}) * h) + omega.a);
    }
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
    const PDEParams &pdep, const SchwarzParams &sp, const SolverParams &solver_params, const Real h
) : PDESolver<Line>(pdep, sp, h), iter_diff(0), iter(0), max_iter(solver_params.max_iter), eps(solver_params.eps) {

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
    for (auto i = 0; i < Nnodes; ++i) {
        u_k[i] = (static_cast<Real>(i)*h - omega.a)*slope + dirichlet.u_a;
    }

#pragma omp parallel for
    for (auto i = 0; i < Nsub; ++i) {
        subdomain_solvers[i].factorize();
    }

    iter_diff = eps + 1;
    while (iter < max_iter && iter_diff > eps) {
        advance();
        iter++;
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

void DiscreteSolver<Line>::print_to_file(std::filesystem::path output) {
    if (!(output.has_filename() && output.has_extension() && output.extension().compare(".vtk") == 0))
        throw std::invalid_argument("The output file or path is not formatted correctly");

    std::string filename = output.c_str();

    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error(
            "Unable to open output file: " + filename
        );
    }

    // Legacy VTK (ASCII) that ParaView can open.
    // We use a 3D rectilinear grid with dimensions (Nnodes, 1, 1).
    out << "# vtk DataFile Version 3.0\n";
    out << "1D solution on uniform mesh\n";
    out << "ASCII\n";
    out << "DATASET RECTILINEAR_GRID\n";
    out << "DIMENSIONS " << Nnodes << " 1 1\n";

    out << std::setprecision(17);

    out << "X_COORDINATES " << Nnodes << " double\n";
    for (Index i = 0; i < Nnodes; ++i) {
        const Real x = omega.a + static_cast<Real>(i) * h;
        out << x;
        if (i + 1 < Nnodes) out << ' ';
        if ((i + 1) % 6 == 0) out << '\n';
    }
    out << "\n";

    out << "Y_COORDINATES 1 double\n0\n";
    out << "Z_COORDINATES 1 double\n0\n";

    out << "POINT_DATA " << Nnodes << "\n";
    out << "SCALARS u double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (Index i = 0; i < Nnodes; ++i) {
        out << u_k[i] << "\n";
    }
}

void DiscreteSolver<Line>::print() const {
    auto& out = std::cout;

    // Legacy VTK (ASCII) that ParaView can open.
    // We use a 3D rectilinear grid with dimensions (Nnodes, 1, 1).
    out << "# vtk DataFile Version 3.0\n";
    out << "1D solution on uniform mesh\n";
    out << "ASCII\n";
    out << "DATASET RECTILINEAR_GRID\n";
    out << "DIMENSIONS " << Nnodes << " 1 1\n";

    out << std::setprecision(17);

    out << "X_COORDINATES " << Nnodes << " double\n";
    for (Index i = 0; i < Nnodes; ++i) {
        const Real x = omega.a + static_cast<Real>(i) * h;
        out << x;
        if (i + 1 < Nnodes) out << ' ';
        if ((i + 1) % 6 == 0) out << '\n';
    }
    out << "\n";

    out << "Y_COORDINATES 1 double\n0\n";
    out << "Z_COORDINATES 1 double\n0\n";

    out << "POINT_DATA " << Nnodes << "\n";
    out << "SCALARS u double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (Index i = 0; i < Nnodes; ++i) {
        out << u_k[i] << "\n";
    }
}

void DiscreteSolver<Line>::advance() {

    Real local_diff = 0.0;
    iter_diff = 0;
    #pragma omp parallel for firstprivate(local_diff)
    for (int i = 0; i < Nsub; ++i) {
        subdomain_solvers[i].update_boundary(current_boundary_cond(i));

        Vector u_i_k = subdomain_solvers[i].solve();

        Boundary non_overlap_bnd = get_subdomain_nonoverlapping_boundary(i);
        Boundary overlap_bnd = get_subdomain_overlapping_boundary(i);

        Index first_node_nonoverlap = get_leftmost_node(non_overlap_bnd);
        Index first_node_overlap = get_leftmost_node(overlap_bnd);

        Index last_node_nonoverlap   = get_rightmost_node(non_overlap_bnd);

        Index local_offset = first_node_nonoverlap - first_node_overlap;

        for (Index k = first_node_nonoverlap; k <= last_node_nonoverlap; ++k) {
            Index local_node = local_offset + (k - first_node_nonoverlap);
            u_next[k] = u_i_k[local_node];
            local_diff += (u_k[k] - u_next[k])*(u_k[k] - u_next[k]);
        }

        local_diff = sqrt(local_diff);
        #pragma omp critical
        if (local_diff >= iter_diff) iter_diff = local_diff;

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
