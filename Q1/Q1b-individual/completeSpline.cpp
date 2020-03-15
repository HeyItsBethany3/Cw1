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

  // Create x values (interpolating points)
  double* x;
  x = new double[n+1];
  const double h = len/double(n);
  for (int i=0; i<n+1; i++) {
    x[i] = i*h;
  }

  // Create f(x) using function pointers
  for (int i=1; i<n; i++) {
    fvec[i]=(*funct)(x[i]);
  }
  fvec[0]=(*funct)(x[0])+((double(1)/double(3))*h*(*derivative)(x[0]));
  fvec[n]=(*funct)(x[n])-((double(1)/double(3))*h*(*derivative)(x[n]));

  // Find diagonal elements
  for (int i=0; i<=n; i++) {
    d[i] = 4;
  }

  // Construct upper diagonal elements
  u[0]=2;
  for (int i=1; i<n; i++) {
    u[i]=1;
  }

  // Construct lower diagonal elements
  l[n-1]=2;
  for (int i=0; i<n-1; i++) {
    l[i]=1;
  }

  delete[] x;
}

// Function prototypes
void findSpline(const double len, const int n, double (*funct)(const double x),
double (*derivative)(const double x), double *d, double *u, double *l, double *fvec);
double f(const double x);
double fD(const double x);


int main(int argc, char* argv[]) {

  double *fvec, *d, *u, *l;
  int n=8;
  fvec = new double[n+1];
  d = new double[n+1];
  u = new double[n];
  l = new double[n];
  // Finds system of equations for spline
  findSpline(2, n, f, fD, d, u, l, fvec);

  // Output d, u, l
  std::cout << "\nd: ";
  for (int i=0; i<n+1; i++) {
    std::cout << d[i] << " ";
  }
  std::cout << "\nu: ";
  for (int i=0; i<n; i++) {
    std::cout << u[i] << " ";
  }
  std::cout << "\nl: ";
  for (int i=0; i<n; i++) {
    std::cout << l[i] << " ";
  }

  // Deallocates storage
  delete[] fvec;
  delete[] d;
  delete[] u;
  delete[] l;

  return 0;
}
