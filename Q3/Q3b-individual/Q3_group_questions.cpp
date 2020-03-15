#include <iostream>
#include <cmath>
#include <fstream> // Used to create file of data
#include <string> // for filename

//FUNCTION PROTOTYPE
double* euler(const double theta0, const double alpha, const double T, const int n, double* y);
double* newtonMethod(const double* x_initial, const double alpha, const double T, const double h, const int &n,
  double* (*InvJ_function)(const double* y_new, const double alpha, const double h),
  double* (*F_function)(const double* y_new, const double* y_old, const double alpha, const double h),
  double* newton_differences, int &k_converge);
void createFiles(const int n, const double* newton_differences, const int k_converge);
double calculateError(double* x_new, double* x_old);
double* invJacobi1(const double* y_new, const double alpha, const double h);
double* func1(const double* y_new, const double* y_old, const double alpha, const double h);



//Beth's Euler

/* Implements euler's method
theta0: initial theta value, T is the time step of interest, alpha is a parameter
of the problem and n is the number of time intervals, y is an empty vector of length 2
Finds vector of (y1,y2) approximate values at time T */
double* euler(const double theta0, const double alpha, const double T, const int n)
{

  const double h = T/double(n); // step-size

  double* newton_differences;
  newton_differences = new double[10]; //creates array of size kmax
  int k_converge;

  double *y_old, *y_new;
  y_old = new double[2];
  y_new = new double[2];

  y_old[0] = theta0;
  y_old[1] = 0;

  for (int i=1; i<=n; i++)
  {
    y_new = newtonMethod(y_old, alpha, T, h, i, invJacobi1, func1, newton_differences, k_converge);

    if (i==1)
    {
      createFiles(i, newton_differences, k_converge);
    }

    y_old[0] = y_new[0];
    y_old[1] = y_new[1];
  }

  delete[] y_old;
  delete[] newton_differences;

  return y_new;
}



// Britta's Newton

double* newtonMethod(const double* x_initial, const double alpha, const double T, const double h, const int &n,
  double* (*InvJ_function)(const double* y_new, const double alpha, const double h),
  double* (*F_function)(const double* y_new, const double* y_old, const double alpha, const double h),
  double* newton_differences, int &k_converge)
{
  double *x_new, *x_old, *y_old;
  x_new = new double[2];
  x_old = new double[2]; // x_old is used as the initial guess for y_new (aka y_(n+1))
  y_old = new double[2]; // y_old refers to y_n, i.e. the result of previous euler step

  x_old[0] = y_old[0] = x_initial[0];
  x_old[1] = y_old[1] = x_initial[1];

  double *Fx, *Inv_Jx;
  Fx = new double[2];
  Inv_Jx = new double[4];

  double error = 1.0; // Arbitrary initialisation
  const double TOL = 1e-12;
  int k = 1; // Counter for number of Netwon iterations
  int kmax = 10;

  while (error > TOL && k<=kmax)
  {
    Fx = (*F_function)(x_old, y_old, alpha, h);
    Inv_Jx = (*InvJ_function)(x_old, alpha, h);

    x_new[0] = x_old[0] - (Inv_Jx[0]*Fx[0] + Inv_Jx[1]*Fx[1]);
    x_new[1] = x_old[1] - (Inv_Jx[2]*Fx[0] + Inv_Jx[3]*Fx[1]);

    error = calculateError(x_new, x_old);

    if (n==1)
    {
      newton_differences[k-1] = calculateError(x_new, x_old);
    }

    x_old[0] = x_new[0];
    x_old[1] = x_new[1];

    //TODO REMOVE PRINT
    //std::cout << "x_" << k << ": (" << x_new[0] << "," << x_new[1] << ")\n";

    k++;
  }

  k_converge = k-1; //The number of iterations it took Newton method to converge

  delete[] Fx;
  delete[] Inv_Jx;
  delete[] x_old;
  delete[] y_old;

  return x_new; //x_new is the final approx for y_new (aka y_(n+1))
}


void createFiles(const int n, const double* newton_differences, const int k_converge)
{
  std::string filename;
  filename = "Newton_n=" + std::to_string(n) + ".csv";
  std::ofstream myfile;
  myfile.open(filename);
  myfile << "k, difference\n";
  for (int j=0; j<k_converge; j++)
  {
    std::cout << "newton_diff[" << j << "] = " << newton_differences[j] << "\n";
    myfile << j+1 << "," << newton_differences[j] << "\n";
  }
  myfile.close();

    //if (n==1 || n == int(T/(4.0*h)) || n == int(T/(2.0*h)))
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



// Q3 group questions functions

double* invJacobi1(const double* y_new, const double alpha, const double h)
{
  double* Inv_Jy;
  Inv_Jy = new double[4];

  double det = 1.0 + pow(h*alpha, 2.0)*cos(y_new[0]);

  Inv_Jy[0] = 1.0 / det;
  Inv_Jy[1] = h / det;
  Inv_Jy[2] = -1.0*h*pow(alpha,2.0)*cos(y_new[0]) / det;
  Inv_Jy[3] = 1.0 / det;

  return Inv_Jy;
}

double* func1(const double* y_new, const double* y_old, const double alpha, const double h)
{
  double* Fy;
  Fy = new double[2];

  Fy[0] = y_new[0] - y_old[0] - h*y_new[1];
  Fy[1] = y_new[1] - y_old[1] + h*pow(alpha,2.0)*sin(y_new[0]);

  return Fy;
}


int main(int argc, char* argv[])
{
  double *y_n;
  y_n = new double[2];

  const double theta0 = M_PI/double(2.0); // initial theta value
  const double alpha = 2.0;
  const double T = 8.0;
  const double h = T/ 32.0;

  const int n = int(T/h);

  y_n = euler(theta0, alpha, T, n);
  // THINK EULER WILL OUTPUT A VECTOR OF Y_N's NOT JUST THE LAST ROOT

  std::cout << "y_n found: (" << y_n[0] << "," << y_n[1] << ")\n";

  return 0;
}
