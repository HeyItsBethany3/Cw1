#include <iostream>
#include <cmath>
#include <cassert>

// FUNCTION PROTOTYPE
double* newtonMethod(const double* x_initial,
  double* (*InvJ_function)(const double* x), double* (*F_function)(const double* x), double &Cmax);
double calculateError(double* x_new, double* x_old);
double verifyConvergence(double* x_new, double* x_old);
double* testF(const double* x);
double* testInvJacobi(const double* x);



double* newtonMethod(const double* x_initial,
  double* (*InvJ_function)(const double* x), double* (*F_function)(const double* x), double &Cmax)
/*
Implements Newton's method
Parameters:
  x_initial: initial guess for the root
  InvJ_function pointer: calculates a vector of entries for the inverse jacobian evaluated at the given x vector
  F_function point: calculates the Newton F(x) function evaluated at given x vector
Returns:
  The root found from Newton's method
*/
{
  double *x_new, *x_old;
  x_new = new double[2];
  x_old = new double[2];

  x_old[0] = x_initial[0];
  x_old[1] = x_initial[1];

  double *Fx, *Inv_Jx;
  Fx = new double[2];
  Inv_Jx = new double[4];

  double C; // quadratic convergence constant for this iteration

  double error = 1.0; // arbitrary initialisation
  const double TOL = 1e-12;
  int k = 1; // counter for number of Netwon iterations
  int kmax = 10; // max number of Newton iterations

  while (error > TOL && k <= kmax)
  {
    Fx = (*F_function)(x_old);
    Inv_Jx = (*InvJ_function)(x_old);

    x_new[0] = x_old[0] - (Inv_Jx[0]*Fx[0] + Inv_Jx[1]*Fx[1]);
    x_new[1] = x_old[1] - (Inv_Jx[2]*Fx[0] + Inv_Jx[3]*Fx[1]);

    error = calculateError(x_new, x_old);

    C = verifyConvergence(x_new, x_old);
    // Finds quadratic convergence constant for this whole example
    if ( C > Cmax)
    {
      Cmax = C;
    }

    x_old[0] = x_new[0];
    x_old[1] = x_new[1];

    k++;
  }

  delete[] Fx;
  delete[] Inv_Jx;
  delete[] x_old;

  return x_new;
}


double calculateError(double* x_new, double* x_old)
/*
Calculates and returns the norm of the difference between the given x_new and x_old vectors.
*/
{
  double error;
  double* error_vector;
  error_vector = new double[2];

  error_vector[0] = x_new[0] - x_old[0];
  error_vector[1] = x_new[1] - x_old[1];

  // Error calculated using vector 2-norm
  error = sqrt(pow(error_vector[0], 2.0) + pow(error_vector[1], 2.0));

  delete[] error_vector;

  return error;
}

double verifyConvergence(double* x_new, double* x_old)
/*
Checks for quadratic convergence to root (1,1), using Definition 6.1
This is done by analysing the convergence constant C = |x_new - root| / |x_old - root|^2
Returns convergence constant
*/
{

  double C; // convergence constant

  double* root; // expected root
  root = new double[2];
  root[0] = root[1] = 1.0;

  C = calculateError(x_new, root) / pow(calculateError(x_old, root), 2.0);

  assert((0 <= C) && (C <= 1));
  //Note: since have exact convergence to the root here needed to allow C=0

  return C;
}


double* testF(const double* x)
/*
Returns F(x) for Lecture 3, slide 8 example.
*/
{
  double *F_output;
  F_output = new double[2];

  F_output[0] = pow(x[0],2.0) - x[1];
  F_output[1] = pow(x[0],2.0) + pow(x[1],2.0) - 2.0;

  return F_output;
}


double* testInvJacobi(const double* x)
/*
Returns Inverse Jacobi evaluated at x for Lecture 3, slide 8 example.
*/
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


int main(int argc, char* argv[])
{
  std::cout.precision(20);

  double *root;
  root = new double[2];
  root[0] = root[1] = 2.0; // Initial guess (2,2)
  double Cmax; // Constant of quadratic convergence

  root = newtonMethod(root, testInvJacobi, testF, Cmax);

  std::cout << "\nAs expected, when we use initial guess (2,2), we find a root = ("
            << root[0] << "," << root[1] << ")\n";
  std::cout << "This had quadratic convergence with constant C = " << Cmax << "\n";

  return 0;
}
