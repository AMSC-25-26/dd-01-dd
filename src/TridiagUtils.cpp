#include <TridiagUtils.hpp>
#include <cmath>

using Real = Types<Line>::Real;
using Index = Types<Line>::Index;
using Size = Types<Line>::Size;
using Vector = Types<Line>::Vector;

FactorizedTridiag::FactorizedTridiag(Size n)
    : _lower(n > 0 ? n - 1 : 0, 0.0),
      _diag(n, 0.0),
      _upper(n > 0 ? n - 1 : 0, 0.0),
      _n(n)
{}

Real& FactorizedTridiag::operator()(Index i, Index j) {
    if (i < 0 || i >= (Index)_n || j < 0 || j >= (Index)_n) {
        throw std::out_of_range("Index out of matrix bounds");
    }

    _is_factorized = false;

    if (i == j) {
        return _diag[i];
    }
    if (j == i + 1) {
        return _upper[i];
    }
    if (j == i - 1) {
        return _lower[j];
    }

    throw std::out_of_range("Index outside tridiagonal band");
}

Real FactorizedTridiag::operator()(Index i, Index j) const {
    if (i < 0 || i >= (Index)_n || j < 0 || j >= (Index)_n) return 0.0;

    if (i == j)     return _diag[i];
    if (j == i + 1) return _upper[i];
    if (j == i - 1) return _lower[j];

    return 0.0;
}

void FactorizedTridiag::factorize() {
    if (_n <= 1) {
        _is_factorized = true;
        return;
    }

    // In-place Thomas algorithm:
    // _diag becomes U diagonal (ûᵢ)
    // _lower becomes L multipliers (mᵢ)

    if (std::abs(_diag[0]) < 1e-15) {
        throw std::runtime_error("Zero pivot encountered at index 0");
    }

    for (Size i = 0; i < _n - 1; ++i) {
        _lower[i] = _lower[i] / _diag[i];                    // mᵢ = lᵢ / ûᵢ
        _diag[i + 1] = _diag[i + 1] - _lower[i] * _upper[i]; // ûᵢ₊₁ = dᵢ₊₁ - mᵢ·uᵢ

        if (std::abs(_diag[i + 1]) < 1e-15) {
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

    Vector y(_n);
    Vector x(_n);

    // Forward substitution: L·y = b
    y[0] = b[0];
    for (Size i = 1; i < _n; ++i) {
        y[i] = b[i] - _lower[i - 1] * y[i - 1];  // _lower now contains mᵢ
    }

    // Backward substitution: U·x = y
    x[_n - 1] = y[_n - 1] / _diag[_n - 1];       // _diag now contains ûᵢ
    for (Index i = (Index)_n - 2; i >= 0; --i) {
        x[i] = (y[i] - _upper[i] * x[i + 1]) / _diag[i];
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
