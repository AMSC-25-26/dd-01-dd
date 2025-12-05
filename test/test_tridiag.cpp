#include <TridiagUtils.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

using Real = Types<1>::Real;
using Index = Types<1>::Index;
using Size = Types<1>::Size;
using Vector = Types<1>::Vector;

Vector multiply(const FactorizedTridiag& A, const Vector& x, Size n) {
    Vector res(n, 0.0);
    for (Index i = 0; i < (Index)n; ++i) {
        Real val = A(i, i) * x[i];
        if (i > 0) val += A(i, i - 1) * x[i - 1];
        if (i < (Index)n - 1) val += A(i, i + 1) * x[i + 1];
        res[i] = val;
    }
    return res;
}

Real calculate_error(const Vector& x, const Vector& x_ref) {
    Real num = 0.0;
    Real den = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        num += (x[i] - x_ref[i]) * (x[i] - x_ref[i]);
        den += x_ref[i] * x_ref[i];
    }
    return std::sqrt(num) / std::sqrt(den);
}

void run_test(Size n, std::mt19937& gen) {
    std::cout << "------------------------------------------------\n";
    std::cout << "Testing with N = " << n << "..." << std::endl;

    FactorizedTridiag A(n);
    std::uniform_real_distribution<> dist_off(-10.0, -1.0); // Off-diagonals (negative)
    std::uniform_real_distribution<> dist_diag(25.0, 50.0); // Diagonal (large positive)
    std::uniform_real_distribution<> dist_sol(-5.0, 5.0);   // Solution values

    for (Index i = 0; i < (Index)n; ++i) {
        Real d = dist_diag(gen);
        Real l = (i > 0) ? dist_off(gen) : 0.0;
        Real u = (i < (Index)n - 1) ? dist_off(gen) : 0.0;

        A(i, i) = d;
        if (i > 0) A(i, i - 1) = l;
        if (i < (Index)n - 1) A(i, i + 1) = u;
    }

    Vector x_exact(n);
    for (Size i = 0; i < n; ++i) x_exact[i] = dist_sol(gen);

    Vector b = multiply(A, x_exact, n);

    try {
        A.factorize();
    } catch (const std::exception& e) {
        std::cerr << "Factorization FAILED: " << e.what() << std::endl;
        return;
    }

    Vector x_calc = A.solve(b);

    Real rel_err = calculate_error(x_calc, x_exact);

    std::cout << "Relative Error: " << std::scientific << rel_err << std::defaultfloat << std::endl;

    if (rel_err < 1e-12) {
        std::cout << "[ OK ] Solution is accurate." << std::endl;
    } else {
        std::cout << "[ FAIL ] Error is too high!" << std::endl;
        exit(1);
    }

    std::cout << "Testing reuse of factorization..." << std::endl;
    for (auto& val : x_exact) val = dist_sol(gen); // New random solution
    Vector b_new = multiply(A, x_exact, n);        // New RHS

    Vector x_calc_new = A.solve(b_new);            // Solve again (cheaply)
    Real rel_err_2 = calculate_error(x_calc_new, x_exact);

    if (rel_err_2 < 1e-12) {
        std::cout << "[ OK ] Reuse successful. Error: " << std::scientific << rel_err_2 << std::endl;
    } else {
        std::cout << "[ FAIL ] Reuse failed!" << std::endl;
        exit(1);
    }
}

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());

    run_test(5, gen);

    run_test(100, gen);

    run_test(1000, gen);

    std::cout << "\nALL TESTS PASSED." << std::endl;
    return 0;
}
