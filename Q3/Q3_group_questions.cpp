#include <iostream>
#include <cmath>
#include <fstream> // Used to create file of data
#include <string> // for filename

//FUNCTION PROTOTYPE
double* euler(const double theta0, const double alpha, const double T, const int n, double* y);
double* newtonMethod(const double* x_initial, const double alpha, const double T, const double h, const int euler_n,
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

  double *newton_differences; //norms of differences for Question 3d
  newton_differences = new double[10]; //creates array of size kmax=10
  int k_converge; //records the number of iterations it takes Newton method to converge

  double *y_old, *y_new;
  y_old = new double[2];
  y_new = new double[2];

  y_old[0] = theta0;
  y_old[1] = 0;

  for (int i=1; i<=n; i++)
  {
    y_new = newtonMethod(y_old, alpha, T, h, i, invJacobi1, func1, newton_differences, k_converge);

    // Outputs for Question 3e
    if (i==1 || i == int(T/(4.0*h)) || i == int(T/(2.0*h)))
    {
      for (int j=0; j<k_converge; j++)
      {
        std::cout << "  For n = " << i << ", newton_diff[" << j << "] = " << newton_differences[j] << "\n";
      }

      // File outputs for Question 3d
      if (n==32.0)
      {
        createFiles(i, newton_differences, k_converge);
      }
    }

    y_old[0] = y_new[0];
    y_old[1] = y_new[1];
  }

  delete[] y_old;
  delete[] newton_differences;

  return y_new;
}



// Britta's Newton

double* newtonMethod(const double* x_initial, const double alpha, const double T, const double h, const int euler_n,
  double* (*InvJ_function)(const double* y_new, const double alpha, const double h),
  double* (*F_function)(const double* y_new, const double* y_old, const double alpha, const double h),
  double* newton_differences, int &k_converge)
  /*
  Implements Newton's method
  Parameters:
    x_initial: initial guess for the root
    alpha: constant to be passed into InvJ and F functions
    h: constant to be passed into InvJ and F functions
    InvJ_function pointer: calculates a vector of entries for the inverse jacobian evaluated at the given y vector
    F_function point: calculates the Newton F(x) function evaluated at given y vector
    euler_n: records which euler iteration we are on - used for outputs for Question 3e
    T: constant, used to identify which euler_n values we want to have outputs for in Question 3e
    newton_differences: vector to hold the norm of the differences between Netwon iterations for Question 3d
    k_converge: keeps track of how many iteratations
  Returns:
    The root found from Newton's method
  */
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
  int kmax = 10; // Max Newton iterations

  while (error > TOL && k<=kmax)
  {
    Fx = (*F_function)(x_old, y_old, alpha, h);
    Inv_Jx = (*InvJ_function)(x_old, alpha, h);

    x_new[0] = x_old[0] - (Inv_Jx[0]*Fx[0] + Inv_Jx[1]*Fx[1]);
    x_new[1] = x_old[1] - (Inv_Jx[2]*Fx[0] + Inv_Jx[3]*Fx[1]);

    error = calculateError(x_new, x_old);

    // Calculates the norm of differences between newton iterations for Question 3d
    if (euler_n==1 || euler_n == int(T/(4.0*h)) || euler_n == int(T/(2.0*h)))
    {
      newton_differences[k-1] = calculateError(x_new, x_old);
    }

    x_old[0] = x_new[0];
    x_old[1] = x_new[1];

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
/*
Creates output csv files for Question 3d.
*/
{
  std::string filename;
  filename = "Q3d_Britta_Newton_n=" + std::to_string(n) + ".csv";
  std::ofstream myfile;
  myfile.open(filename);
  assert(myfile.is_open());
  myfile << "k, difference\n";
  for (int j=0; j<k_converge; j++)
  {
    myfile << j+1 << "," << newton_differences[j] << "\n";
  }
  myfile.close();
  std::string command = "mv "+filename+" Documents/GitHub/Cw1/Q3/";
  system(command.c_str());
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



// Q3 group questions functions - Beth and Britta's functions

double* invJacobi1(const double* y_new, const double alpha, const double h)
/*
Returns Inverse Jacobi evaluated at x for Question 3.
*/
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
/*
Returns F(x) for Question 3.
*/
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

  //Question 3e
  double h = T; //initialise h
  int n;
  for (int i=1; i<=20; i++)
  {
    h = h/2.0; //Note: this will also check for h=T/32.0 which is needed for Question 3d
    n = int(T/h);
    std::cout << "\n\nFor h = " << T/double(n) << "\n";
    y_n = euler(theta0, alpha, T, n);
  }

  return 0;
}
