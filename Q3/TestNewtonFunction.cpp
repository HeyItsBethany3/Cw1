#include "TestNewtonFunction.hpp"
#include "TestJacobiInvMatrix.hpp"
#include <cmath>

// Constructor
TestNewtonFunction::TestNewtonFunction()
{
}

// Destructor
TestNewtonFunction::~TestNewtonFunction()
{
}

// Overide pure virtual method
TestNewtonFunction::double* evaluateF(const double* x)
/*
F(x) from Lecture 3, slide 8.
*/
{
  double* F_output;
  F_output = new double[2];

  F_output[0] = pow(x[0],2.0) - x[1];
  F_output[1] = pow(x[0],2.0) + pow(x[1],2.0) - 2.0;

  return F_output;
}

// Overide pure virtual method
TestNewtonFunction::AbstractMatrix evaluateInvJacobi(const double* x)
{
  TestJacobiInvMatrix J_output(x);
  return J_output;
}
