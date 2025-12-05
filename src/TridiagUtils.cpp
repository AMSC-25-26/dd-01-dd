#include <TridiagUtils.hpp>
#include <cmath> // for std::abs

using Real = Types<1>::Real;
using Index = Types<1>::Index;
using Size = Types<1>::Size;
using Vector = Types<1>::Vector;

FactorizedTridiag::FactorizedTridiag(Size n)
    : _lower(n > 0 ? n - 1 : 0, 0.0),
      _diag(n, 0.0),
      _upper(n > 0 ? n - 1 : 0, 0.0),
      _n(n),
      _du(n, 0.0),
      _dl(n > 0 ? n - 1 : 0, 0.0)
{}

// Read-write access
Real& FactorizedTridiag::operator()(Index i, Index j) {
    if (i < 0 || i >= (Index)_n || j < 0 || j >= (Index)_n) {
        throw std::out_of_range("Index out of matrix bounds");
    }

    // IMPORTANT: If we modify the matrix, the previous factorization is wrong.
    _is_factorized = false;

    if (i == j) {
        return _diag[i];
    }
    if (j == i + 1) {
        return _upper[i];
    }
    if (j == i - 1) {
        return _lower[j]
    }

    throw std::out_of_range("Index outside tridiagonal band");
}

// Read-only access
Real FactorizedTridiag::operator()(Index i, Index j) const {
    if (i < 0 || i >= (Index)_n || j < 0 || j >= (Index)_n) return 0.0;

    if (i == j)     return _diag[i];
    if (j == i + 1) return _upper[i];
    if (j == i - 1) return _lower[j];

    return 0.0; // Elements outside the band are zero
}

void FactorizedTridiag::factorize() {
    if (_n == 0) return;

    // Thomas Algorithm (LU Decomposition)
    // We compute:
    //   _du[i] (diagonal of U)
    //   _dl[i] (lower diagonal of L)
    // The upper diagonal of U is the same as the original _upper.

    // 1. First element
    _du[0] = _diag[0];
    if (std::abs(_du[0]) < 1e-15) {
        throw std::runtime_error("Zero pivot encountered at index 0");
    }

    // 2. Loop for remaining elements
    for (Size i = 0; i < _n - 1; ++i) {
        // Calculate L multiplier
        _dl[i] = _lower[i] / _du[i];

        // Calculate new U diagonal
        _du[i + 1] = _diag[i + 1] - _dl[i] * _upper[i];

        // Stability check
        if (std::abs(_du[i + 1]) < 1e-15) {
            throw std::runtime_error("Zero pivot encountered during factorization");
        }
    }

    _is_factorized = true;
}

Vector FactorizedTridiag::solve(const Vector& b) const {
    if (!_is_factorized) {
        throw std::logic_error("Matrix is not factorized. Call factorize() first.");
    }
    if (b.size() != _n) {
        throw std::invalid_argument("Vector b size must match matrix size");
    }

    Vector y(_n); // Intermediate vector
    Vector x(_n); // Solution vector

    // 1. Forward substitution: L * y = b
    // L is unit lower triangular: diag is 1, sub-diag is _dl
    y[0] = b[0];
    for (Size i = 1; i < _n; ++i) {
        y[i] = b[i] - _dl[i - 1] * y[i - 1];
    }

    // 2. Backward substitution: U * x = y
    // U has diagonal _du and super-diagonal _upper
    x[_n - 1] = y[_n - 1] / _du[_n - 1];

    for (Index i = (Index)_n - 2; i >= 0; --i) {
        x[i] = (y[i] - _upper[i] * x[i + 1]) / _du[i];
    }

    return x;
}

void FactorizedTridiag::print() const {
    std::cout << "Tridiagonal Matrix (" << _n << "x" << _n << "):\n";
    for (Size i = 0; i < _n; i++) {
        for (Size j = 0; j < _n; j++) {
            std::cout << this->operator()(i, j) << "\t";
        }
        std::cout << "\n";
    }
}
