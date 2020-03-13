#include <iostream>
#include <cmath>


double* newtonMethod(const double* x, double* (*J_function)(const double* x), double* (*F_function)(const double* x))
{
  double *x_new, *x_old;
  x_new = new double[2];
  x_old = new double[2];

  x_old = x; //TODO check if this assignment works

  double error = 1;
  const double TOL = 1e-12;

  double F_1, F_2;
  //TODO MAKE THE SAME FOR J

  while (error > TOL)
  {

  }


}

double calculateError(double* x_new, double* x_old)
{
  double error;
  double* error_vector;
  error_vector = new double[2];

  error_vector[0] = x_new[0] - x_old[0];
  error_vector[1] = x_new[1] - x_old[1];

  //TODO used Euclidean Norm - is this right? or meant to use L2-norm?
  error = sqrt(pow(error_vector[0], 2.0) + pow(error_vector[1], 2.0));

  delete[] error_vector;

  return error;
}


double* testF(const double* x)
{
  double *F_output;
  F_output = new double[2];

  F_output[0] = pow(x[0],2.0) - x[1];
  F_output[1] = pow(x[0],2.0) + pow(x[1],2.0) - 2.0;

  return F_output;
}

double* testInvJacobi(const double* x)
{
  double *J_output;
  J_output = new double[4];

  double det;

  det = 4.0*x[0]*x[1] + 2.0*x[0];

  J_output[0] = 2.0*x[1] / det;
  J_output[1] = 1.0 / det;
  J_output[2] = -2.0*x[0] / det;
  J_output[3] = 2.0*x[0] / det;

  return J_output;
}
