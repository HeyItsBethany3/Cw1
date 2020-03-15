#include <iostream>
#include <cmath>
using namespace std;

double Gaussian_3pnt(int n);
double Gaussian_5pnt(int n);
double f(double x, int n);

double Gaussian_3pnt(int n)
{
    double x0, x1, x2;

    x0 = sqrt(3) /sqrt(5);
    x1 = 0.;
    x2 = - x0;

    double sum = (8./9.)*f(x1,n) + (5./9.)*f(x0,n) + (5./9.)*f(x2,n);
    return sum;
}

double Gaussian_5pnt(int n)
{
    double x0, x1, x2, x3, x4;

    x0 = 0.;
    x1 = (1./3.) * sqrt(5. - 2. * sqrt(10./7.));
    x2 = (-1./3.) * sqrt(5. - 2. * sqrt(10./7.));
    x3 = (1./3.) * sqrt(5. + 2.*sqrt(10./7.));
    x4 = (-1./3.) * sqrt(5. + 2.*sqrt(10./7.));

    double sum = (128./225.) * f(x0,n)
               + ((322. + 13. * sqrt(70.))/900.) * (f(x1,n) + f(x2,n))
               + ((322. - 13. * sqrt(70.))/900.) * (f(x3,n) + f(x4,n));

    return sum;
}

// function to test
double f(double x, int n)
{
    return pow(x, n);
}

int main(int argc, char* argv[])
{
    cout << "~~~~~~~~3 point Gaussian quadrature rule~~~~~~~~\n";
    for (int k=0; k<=10; k++)
        {
            cout << "\nFor k = " << k << " : " << Gaussian_3pnt (k);
        }

    cout << "\n\n~~~~~~~~5 point Gaussian quadrature rule~~~~~~~~\n";
    for (int k=0; k<=10; k++)
        {
            cout << "\nFor k = " << k << " : " << Gaussian_5pnt (k);
        }
    return 0;
}
