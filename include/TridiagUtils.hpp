#include <types.hpp>
#include <iostream>

class FactorizedTridiag  : protected Types{
private:
    Vector lower_;  // n-1 elements: l[0] ... l[n-2]
    Vector diag_;   // n elements:   d[0] ... d[n-1]
    Vector upper_;  // n-1 elements: u[0] ... u[n-2]
    Size n_;

public:
    FactorizedTridiag(Size n);

    Real& operator()(Index i, Index j);

    Vector solve(const Vector& b);

    void print() const;
};
