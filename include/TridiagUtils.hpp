#ifndef DD_01_DD_TRIDIAGUTILS_HPP
#define DD_01_DD_TRIDIAGUTILS_HPP

#include <types.hpp>
#include <iostream>
#include <vector>
#include <stdexcept>

class FactorizedTridiag  : protected Types<Line> {
private:
    Vector _lower;  // Lower diagonal (n-1 elements)
    Vector _diag;   // Main diagonal (n elements)
    Vector _upper;  // Upper diagonal (n-1 elements)
    Size _n;        // Dimension of the matrix

    // --- LU Factorization Data ---
    // U is upper triangular with diagonal _du and super-diagonal _upper (original)
    // L is unit lower triangular with sub-diagonal multipliers _dl
    Vector _du;    
    Vector _dl;     

    bool _is_factorized = false; // Flag to ensure solve() is not called before factorize()

public:
    explicit FactorizedTridiag(Size n);

    // Read-write access. modifying the matrix invalidates the factorization.
    Real& operator()(Index i, Index j);

    // Read-only access.
    Real operator()(Index i, Index j) const;

    void factorize();

    Vector solve(const Vector& b) const;

    void print() const;
};

#endif //DD_01_DD_TRIDIAGUTILS_HPP
