#include <iostream>
#include <cmath>


// f(x)
double f(const double x) {
  return (exp(x)* cos(1.25*M_PI*x));
}


// f'(x)
double fD(const double x) {
  return (exp(x)*(cos(1.25*M_PI*x)-(1.25*M_PI*sin(1.25*M_PI*x))));
}


/* Function which creates a complete spline interpolant for a function f on a
grid between 0 and len with step-size h=len/n (there are n+1 interpolating points)
It finds the tridiagonal system of equations which can be solved to find the
coefficients of the spline
It returns vectors d (diagonal), u (upper diagonal), l (lower diagonal),
f (fvalues) and n+1 (length of vectors)
*/
void findSpline(const double len, const int n, double (*funct)(const double x),
double (*derivative)(const double x), double *d, double *u, double *l, double *fvec) {
  //TODO: Do we approximate the derivative f' or use exact?

  // Create x values (interpolating points)
  double* x;
  x = new double[n+1];
  const double h = len/double(n);
  for (int i=0; i<n+1; i++) {
    x[i] = i*h;
  }

  // Create f(x) using function pointers
  fvec = new double[n+1];
  for (int i=1; i<n; i++) {
    fvec[i]=(*funct)(x[i]);
  }
  fvec[0]=(*funct)(x[0])+((double(1)/double(3))*h*(*derivative)(x[0]));
  fvec[n]=(*funct)(x[n])-((double(1)/double(3))*h*(*derivative)(x[n]));



  // Find diagonal elements
  d = new double[n+1];
  for (int i=0; i<=n; i++) {
    d[i] = 4;
  }

  // Construct upper diagonal elements
  u = new double[n+1];
  u[0]=2;
  for (int i=1; i<=n; i++) {
    u[i]=1;
  }

  // Construct lower diagonal elements
  l = new double[n+1];
  l[n]=2;
  for (int i=0; i<n; i++) {
    l[i]=1;
  }

  delete[] x;

}


int main(int argc, char* argv[]) {

  double *fvec, *d, *u, *l;
  findSpline(10, 5, f, fD, d, u, l, fvec);


  delete[] fvec;
  delete[] d;
  delete[] u;
  delete[] l;

  return 0;
}
// TODO: create header file etc
