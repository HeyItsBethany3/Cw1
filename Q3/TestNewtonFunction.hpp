#ifndef TESTNEWTONFUNCTION
#define TESTNEWTONFUNCTION

#include "AbstractNewtonFunction.hpp"

class TestNewtonFunction: public AbstractNewtonFunction
/*
System from Lecture 3, slide 8.
*/
{
  public:
    // Constructor
    TestNewtonFunction();

    // Destructor
    ~TestNewtonFunction();

    // Overide pure virtual method
    double* evaluateF(const double* x);

    // Overide pure virtual method
    AbstractMatrix evaluateInvJacobi(const double* x);

}


#endif
