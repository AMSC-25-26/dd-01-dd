#include <TridiagUtils.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

using Real = Types<Line>::Real;
using Index = Types<Line>::Index;
using Size = Types<Line>::Size;
using Vector = Types<Line>::Vector;

// Multiply tridiagonal matrix by vector: res = A * x
// Must be called BEFORE factorize() since factorize() modifies the matrix
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
    std::uniform_real_distribution<> dist_off(-10.0, -1.0);
    std::uniform_real_distribution<> dist_diag(25.0, 50.0);
    std::uniform_real_distribution<> dist_sol(-5.0, 5.0);

    // Fill matrix
    for (Index i = 0; i < (Index)n; ++i) {
        Real d = dist_diag(gen);
        Real l = (i > 0) ? dist_off(gen) : 0.0;
        Real u = (i < (Index)n - 1) ? dist_off(gen) : 0.0;

        A(i, i) = d;
        if (i > 0) A(i, i - 1) = l;
        if (i < (Index)n - 1) A(i, i + 1) = u;
    }

    // Generate exact solution and compute RHS BEFORE factorization
    Vector x_exact(n);
    for (Size i = 0; i < n; ++i) x_exact[i] = dist_sol(gen);
    Vector b = multiply(A, x_exact, n);

    // Also prepare second test case BEFORE factorization
    Vector x_exact_2(n);
    for (Size i = 0; i < n; ++i) x_exact_2[i] = dist_sol(gen);
    Vector b_2 = multiply(A, x_exact_2, n);

    // Now factorize (modifies A in-place)
    try {
        A.factorize();
    } catch (const std::exception& e) {
        std::cerr << "Factorization FAILED: " << e.what() << std::endl;
        return;
    }

    // Test 1: solve first system
    Vector x_calc = A.solve(b);
    Real rel_err = calculate_error(x_calc, x_exact);

    std::cout << "Relative Error: " << std::scientific << rel_err << std::defaultfloat << std::endl;

    if (rel_err < 1e-12) {
        std::cout << "[ OK ] Solution is accurate." << std::endl;
    } else {
        std::cout << "[ FAIL ] Error is too high!" << std::endl;
        exit(1);
    }

    // Test 2: reuse factorization with different RHS
    std::cout << "Testing reuse of factorization..." << std::endl;
    Vector x_calc_2 = A.solve(b_2);
    Real rel_err_2 = calculate_error(x_calc_2, x_exact_2);

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
