#include <iostream>
#include <cmath>

using namespace std;

double* newton(void(*F)(double, double, double*), void(*dF)(double, double, double**), double* x0, double TOL);
double** allocate_matrix();
void deallocate_matrix(double** matrix);
double* matrixxvector(double** mat, double* vec);
void matrixinvert(double** matrix, double** invert);
double norm(double a, double b);
void F(double x1, double x2, double* y);
void dF(double x1, double x2, double** y);
double convergence(double* x, double* new_x, double* root, int order);

// Newton function
double* newton(void(*F)(double, double, double*), void(*dF)(double, double, double**), double* x0, double TOL)
{
    int counter = 0;
    int nmax = 100;
    double diff;
    double conv;

    double* new_x;
    new_x =  new double[2];

    double* x;
    x = new double[2];

    double* update;
    update =  new double[2];

    double* fun;
    fun = new double[2];

    double** Jacobian;
    Jacobian = allocate_matrix();

    double** InvertJacobian;
    InvertJacobian = allocate_matrix();

    double* root;
    root = new double[2];

    x[0] = x0[0];
    x[1] = x0[1];

    do
        {
            dF(x[0], x[1], Jacobian);
            F(x[0], x[1], fun);

            matrixinvert(Jacobian, InvertJacobian);

            update = matrixxvector(InvertJacobian, fun);

            new_x[0] = x[0] + update[0];
            new_x[1] = x[1] + update[1];

            diff = norm(new_x[0] - x[0], new_x[1] - x[1]);

            conv = convergence(x, new_x, root, 2);
            if(conv < 1 || conv == 1)
            {
                cout << "Newton's order of convergence is indeed quadratic \n";
            }
            else
            {
                cout << " Wrong convergence. \n";
            }

            x[0] = new_x[0];
            x[1] = new_x[1];
            counter++;

            cout << "and the updated root guess is: (" << x[0] << ", " << x[1] << ")\n\n";
        }
    while(diff > TOL && counter <= nmax);

    delete[] x;
    delete[] update;
    delete[] fun;

    deallocate_matrix(Jacobian);
    deallocate_matrix(InvertJacobian);

    return new_x;
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
            cout << "Matrix is not invertible.";
        }
}

// function to find the norm of 2 vectors
double norm(double a, double b)
{
    return pow((pow(a, 2.) + pow(b, 2.)), 0.5);
}

// function to test
void F(double x1, double x2, double* y)
{
    y[0] = pow(x1, 2.) - x2;
    y[1] = pow(x1, 2.) + pow(x2, 2.) - 2;
}

// derivative of function to test
void dF(double x1, double x2, double** y)
{
    y[0][0] = 2 * x1;
    y[0][1] = -1;
    y[1][0] = 2 * x1;
    y[1][1] = 2 * x2;
}

// function to check the convergence
double convergence(double* x, double* new_x, double* root, int order)
{
    double numenator = norm(new_x[0] - root[0], new_x[1] - root[1]);
    double denominator = pow(norm(x[0] - root[0], x[1] - root[1]), order);

    return numenator/denominator;
}

int main(int argc, char* argv[])
{
    int counter;
    double TOL = 1e-10;
    double* root;
    root = new double[2];

    double* x0;
    x0 = new double[2];

    x0[0] = 2.;
    x0[1] = 2.;

    root = newton(F, dF, x0, TOL);
    cout << "Finally, after " << counter << " iterations, Newton's root is: ("<< root[0] << ", "<< root[1] << ") \n";

    return 0;
}
