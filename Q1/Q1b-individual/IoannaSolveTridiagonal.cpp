#include <iostream>
using namespace std;

double* solve_Tridiagonal(const int n, double *l, double *d, double *u, double *f, char &flag);
void test(double *l, double *u, double *d, double *f, char &flag);

// function to give the solution of tridiagonal systems
double* solve_Tridiagonal(const int n, double *l, double *d, double *u, double *f, char &flag)
{
    double* x = new double[n+1];
    flag = 't';

    // elimination step
    for(int i=1; i<=n; i++)
        {
            d[i] -= u[i-1] * (l[i] / d[i-1]);
            f[i] -= f[i-1] * (l[i] / d[i-1]);
            if(d[i] == 0)
            {
                flag = 'f';
                break;
            }
        }

    if(flag == 't')
    {
        // backsolve step
        x[n] = f[n]/d[n];
        for(int i=n-1; i>=0; i--)
            {
                x[i] = (f[i] - u[i] * x[i+1]) / d[i];
            }

    return x;
    }

    else
    {
        cout << "d = 0" << "\n";
    }
}

// function to test with exercise 3
void test(double *l, double *d, double *u, double *f, char &flag)
{
    int n = 3;
    double* solution;

    l = new double[n+1];
    l[0] = 0.0;
    l[1] = 2.0;
    l[2] = 1.0;
    l[3] = 1.0;

    d = new double[n+1];
    d[0] = 6.0;
    d[1] = 4.0;
    d[2] = 4.0;
    d[3] = 6.0;

    u = new double[n+1];
    u[0] = 1.0;
    u[1] = 1.0;
    u[2] = 2.0;
    u[3] = 0.0;

    f = new double[n+1];
    f[0] = 8.0;
    f[1] = 13.0;
    f[2] = 22.0;
    f[3] = 27.0;

    solution = new double[n+1];
    solution = solve_Tridiagonal(n, l, d, u, f, flag);

    if(flag == 't')
    {
        cout << "The solution of the tridiagonal matrix is: " << "\n";
        for(int i=0; i<=n; i++)
            {
                cout << "x[" << i << "] = " << solution[i] << "\n";
            }
    }
    delete[] solution;
}

int main(int argc, char* argv[])
{
    double *l, *d, *u, *f;
    char flag = 't';

    test(l, d, u, f, flag);

    delete[] l;
    delete[] d;
    delete[] u;
    delete[] f;

    return 0;
}
