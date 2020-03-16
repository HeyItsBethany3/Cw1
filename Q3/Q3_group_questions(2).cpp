


#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

// function protypes
double* euler(const double theta0, const double alpha, const double T, const int n);
double* newtonMethod(void(*F)(double, double, double, double, double*, const double, const double),
 void(*dF)(double, double, double**, const double, const double), double* x0,
 const double alpha, const double T, const double h, const int &n,
 double* newton_differences, int &k_converge);
double** allocate_matrix();
void deallocate_matrix(double** matrix);
double* matrixxvector(double** mat, double* vec);
void matrixinvert(double** matrix, double** invert);
double norm(double a, double b);
void F(double y_new1, double y_new2, double y_old1, double y_old2, double* Fy, const double alpha, const double h);
void dF(double y_new1, double y_new2, double** y, const double alpha, const double h);


// Euler's method
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
        y_new = newtonMethod(F, dF, y_old, alpha, T, h, i, newton_differences, k_converge);

        // Outputs for Question 3e
        if (i==1 || i == int(T/(4.0*h)) || i == int(T/(2.0*h)))
        {
            for (int j=0; j<k_converge; j++)
            {
                std::cout << "  For n = " << i << ", newton_diff[" << j << "] = " << newton_differences[j] << "\n";
            }
        }

        y_old[0] = y_new[0];
        y_old[1] = y_new[1];
    }

        delete[] y_old;
        delete[] newton_differences;

        return y_new;
}

// Newton's method
double* newtonMethod(void(*F)(double, double, double, double, double*, const double, const double),
 void(*dF)(double, double, double**, const double, const double), double* x0,
 const double alpha, const double T, const double h, const int &n,
 double* newton_differences, int &k_converge)
{
    int k = 1, nmax = 10;
    double diff = 1., conv;
    double const TOL = 1e-12;

    double *x, * new_x, *update, *Fy, *y;
    x = new double[2];
    new_x = new double[2];
    update =  new double[2];
    Fy = new double[2];
    y = new double[2];

    x[0] = y[0] = x0[0];
    x[1] = y[1] = x0[1];

    double **Jacobian, **InvertJacobian;;
    Jacobian = allocate_matrix();
    InvertJacobian = allocate_matrix();

    while(diff > TOL && k <= nmax)
        {
            dF(x[0], x[1], Jacobian, alpha, h);
            F(x[0], x[1], y[0], y[1], Fy, alpha, h);

            matrixinvert(Jacobian, InvertJacobian);

            update = matrixxvector(InvertJacobian, Fy);

            new_x[0] = x[0] + update[0];
            new_x[1] = x[1] + update[1];

            diff = norm(new_x[0] - x[0], new_x[1] - x[1]);

            if (n==1 || n == int(T/(4.0*h)) || n == int(T/(2.0*h)))
            {
              newton_differences[k] = diff;
            }

            x[0] = new_x[0];
            x[1] = new_x[1];

            k++;
        }

    k_converge = k-1;

    delete[] x;
    delete[] y;
    delete[] update;
    delete[] Fy;

    deallocate_matrix(Jacobian);
    deallocate_matrix(InvertJacobian);

    return new_x;
}

// function to test
void F(double y_new1, double y_new2, double y_old1, double y_old2, double* Fun, const double alpha, const double h)
{
    double* Fy;
    Fy = new double[2];

    Fy[0] = y_new1 - y_old1 - h * y_new2;
    Fy[1] = y_new2 - y_old2 + h * pow(alpha, 2.) * sin(y_new1);
}

// derivative of function to test
void dF(double y_new1, double y_new2, double** y, const double alpha, const double h)
{
    y[0][0] = 1.;
    y[0][1] = -h;
    y[1][0] =  h * pow(alpha, 2.) * cos(y_new1);
    y[1][1] = 1.;
}

// function to allocate storage for 2x2 matrix
double** allocate_matrix()
{
    double** matrix;
    matrix = new double*[2];
    for (int i=0; i<2; i++)
        {
            matrix[i] = new double[2];
        }

    return matrix;
}

// function to free the storage used
void deallocate_matrix(double** matrix)
{
    delete[] matrix[0];
    delete[] matrix[1];
    delete[] matrix;
}

// function to multiply a 2x2 matrix with 2x1 vector
double* matrixxvector(double** mat, double* vec)
{
    double* matvec;
    matvec = new double[2];

    // (-1) is for the -F in Newton
    matvec[0] = -1 * (mat[0][0] * vec[0] + mat[0][1] * vec[1]);
    matvec[1] = -1 * (mat[1][0] * vec[0] + mat[1][1] * vec[1]);

    return matvec;
}

// function to invert a 2x2 matrix
void matrixinvert(double** matrix, double** invert)
{
    double denominator = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    if(denominator != 0)
        {
            invert[0][0] = matrix[1][1] / denominator;
            invert[0][1] = (-1. * matrix[0][1]) / denominator;
            invert[1][0] = (-1. * matrix[1][0]) / denominator;
            invert[1][1] = matrix[0][0] / denominator;
        }
    else
        {
            std::cout << "Matrix is not invertible.";
        }
}

// function to find the norm of 2 vectors
double norm(double a, double b)
{
    return pow((pow(a, 2.) + pow(b, 2.)), 0.5);
}


int main(int argc, char* argv[])
{
  double *y_n;
  y_n = new double[2];

  const double theta0 = M_PI / double(2.);
  const double alpha = 2.0;
  const double T = 8.0;

  //Question 3e
  double h = T;
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


