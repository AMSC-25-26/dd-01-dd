#include "types.hpp"
#include <iostream>

class FactorizedTridiag {
private:
    Types::Vector lower_;  // n-1 elements: l[0] ... l[n-2]
    Types::Vector diag_;   // n elements:   d[0] ... d[n-1]
    Types::Vector upper_;  // n-1 elements: u[0] ... u[n-2]
    Types::Size n_;

public:
    FactorizedTridiag(Types::Size n);

    Types::Real& operator()(Types::Index i, Types::Index j);

    Types::Vector solve(const Types::Vector& b);

    void print() const;
};
