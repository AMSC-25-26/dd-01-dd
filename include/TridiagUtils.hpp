#ifndef DD_01_DD_TRIDIAGUTILS_HPP
#define DD_01_DD_TRIDIAGUTILS_HPP

#include <types.hpp>
#include <iostream>

class FactorizedTridiag  : protected Types<Line> {
private:
    Vector _lower;  // n-1 elements: l[0] ... l[n-2]
    Vector _diag;   // n elements:   d[0] ... d[n-1]
    Vector _upper;  // n-1 elements: u[0] ... u[n-2]
    Size _n;

public:
    FactorizedTridiag(Size n);

    Real& operator()(Index i, Index j);

    Vector solve(const Vector& b);

    void print() const;
};

#endif //DD_01_DD_TRIDIAGUTILS_HPP
