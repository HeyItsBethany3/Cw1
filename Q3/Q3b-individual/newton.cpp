#include <iostream>
#include <cmath>

// FUNCTION PROTOTYPE
double* newtonMethod(const double* x, double* (*InvJ_function)(const double* x), double* (*F_function)(const double* x));
double calculateError(double* x_new, double* x_old);
double* testF(const double* x);
double* testInvJacobi(const double* x);



double* newtonMethod(const double* x, double* (*InvJ_function)(const double* x), double* (*F_function)(const double* x))
{
  double *x_new, *x_old;
  x_new = new double[2];
  x_old = new double[2];

  x_old[0] = x[0];
  x_old[1] = x[1];
  std::cout << "x_0: (" << x_old[0] << " , " << x_old[1] << ")\n";

  double error = 1.0;
  const double TOL = 1e-12;
  double *Fx, *Inv_Jx;
  Fx = new double[2];
  Inv_Jx = new double[4];
  int k=1;
  while (error > TOL)
  {
    Fx = (*F_function)(x_old);
    Inv_Jx = (*InvJ_function)(x_old);

    x_new[0] = x_old[0] - (Inv_Jx[0]*Fx[0] + Inv_Jx[1]*Fx[1]);
    x_new[1] = x_old[1] - (Inv_Jx[2]*Fx[0] + Inv_Jx[3]*Fx[1]);

    error = calculateError(x_new, x_old);
    x_old = x_new;

    std::cout << "x_" << k << ": (" << x_new[0] << " , " << x_new[1] << ")\n";
    std::cout << "  error = " << error << "\n";
    std::cout << "  Fx_" << k-1 << ": (" << Fx[0] << " , " << Fx[1] << ")\n\n";
    //TODO WHY IS THIS Fx DIFFERENT TO WHAT IS CALCULATED IN TESTF FUNCTION?

    k++;
  }


  delete[] Fx;
  delete[] Inv_Jx;
  delete[] x_old;

  return x_new;
}


double calculateError(double* x_new, double* x_old)
{
  double error;
  double* error_vector;
  error_vector = new double[2];

  error_vector[0] = x_new[0] - x_old[0];
  error_vector[1] = x_new[1] - x_old[1];

  // Error calculated using Euclidean norm
  error = sqrt(pow(error_vector[0], 2.0) + pow(error_vector[1], 2.0));

  delete[] error_vector;

  return error;
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

  //TODO REMOVE PRINTS - looks like correctly calculated
  std::cout << "x: (" << x[0] << " , " << x[1] << ")\n";
  std::cout << "Fx: (" << F_output[0] << " , " << F_output[1] << ")\n";

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
  double *root;
  root = new double[2];

  // Intial guess
  root[0] = root[1] = 2.0;

  root = newtonMethod(root, testF, testInvJacobi);

  std::cout << "\nroot = (" << root[0] << " , " << root[1] << ")\n";


  return 0;
}
