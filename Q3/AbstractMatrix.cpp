#include "AbstractMatrix.hpp"

// Destructor
AbstractMatrix::~AbstractMatrix()
{
}

//Matrix-vector multiplication
double* AbstractMatrix::operator*(const double* v) const
{
  double* answer;
  answer = new double[2];

  answer[0] = mData[0][0]*v[0] + mData[0][1]*v[1]
  answer[1] = mData[1][0]*v[0] + mData[1][1]*v[1]

  return answer;
}
