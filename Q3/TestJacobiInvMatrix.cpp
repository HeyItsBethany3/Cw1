#include "TestJacobiInvMatrix.hpp"
#include<cmath>

// Specialised Constructor
TestJacobiInvMatrix::TestJacobiInvMatrix(const double* x)
{
  double det;

  det = 4.0*x[0]*x[1] + 2.0*x[0];

  mData[0][0] = 2.0*x[1] / det;
  mData[0][1] = 1.0 / det;
  mData[1][0] = -2.0*x[0] / det;
  mData[1][1] = 2.0*x[0] / det;
}



// Destructor
TestJacobiInvMatrix::~TestJacobiInvMatrix()
{
}
