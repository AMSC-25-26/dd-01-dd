#ifndef DD_01_DD_TRIDIAGUTILS_HPP
#define DD_01_DD_TRIDIAGUTILS_HPP

#include <types.hpp>
#include <iostream>
#include <vector>
#include <stdexcept>

class FactorizedTridiag  : protected Types<Line> {
private:
    Vector _lower;  
    Vector _diag;  
    Vector _upper;
    Size _n;

    bool _is_factorized = false;

public:
    FactorizedTridiag(Size n);

    // Read-write access. modifying the matrix invalidates the factorization.
    Real& operator()(Index i, Index j);

    // Read-only access.
    Real operator()(Index i, Index j) const;

    void factorize();

    Vector solve(const Vector& b) const;

    void print() const;
};

#endif //DD_01_DD_TRIDIAGUTILS_HPP
