#ifndef ABSTRACTNEWTONFUNCTION
#define ABSTRACTNEWTONFUNCTION

#include "AbstractMatrix.hpp"

class AbstractNewtonFunction {

  public:
    // Destructor
    virtual ~AbstractNewtonFunction() = 0;

    // F(x)
    virtual double* evaluateF(const double* x) = 0;

    // Inverse J(x)
    virtual AbstractMatrix evaluateInvJacobi(const double* x) = 0;

    // Inverse Jacobi
    //AbstractMatrix* mInvJacobi;   //TODO IS THIS NEEDED? already made evaluateInvJacobi return the matrix

};

#endif
