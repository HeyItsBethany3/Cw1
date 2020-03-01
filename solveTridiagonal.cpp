#include <iostream>

//FUNCTION PROTOTYPE
double* solveTridiaognal(const int n, double *d, double *u, double *l, double *fvec);
void testExample(double *d, double *u, double *l, double *fvec);




double* solveTridiaognal(const int n, double *d, double *u, double *l, double *fvec)
/*
Takes a complete spline tridiagonal system as input and solves it.
Returns the B-spline coefficients.
*/
{

  double* c; //B-spline coefficients
  c = new double[n+1];

  // Elimination stage
  for(int i=1; i<=n; i++)
  {
    d[i] = d[i] - u[i-1]*(l[i-1]/d[i-1]);
    fvec[i] = fvec[i] - fvec[i-1]*(l[i-1]/d[i-1]);
  }

  //Backsolve
  c[n] = fvec[n]/d[n];
  for(int i=n-1; i>=0; i--)
  {
    c[i] = ( fvec[i] - u[i]*c[i+1] )/d[i];
  }

  return c;
}



void testExample(double *d, double *u, double *l, double *fvec)
/*
Exercise 3, Section 2.6, Epperson 2013
*/
{
  int n=3;
  double* answer;

  d = new double[n+1];
  d[0] = d[3] = 6.0;
  d[1] = d[2] = 4.0;

  u = new double[n+1];
  u[0] = u[1] = 1.0;
  u[2] = 2.0;


  l = new double[n+1];
  l[0] = 2.0;
  l[1] = l[2] = 1.0;

  fvec = new double[n+1];
  fvec[0] = 8.0;
  fvec[1] = 13.0;
  fvec[2] = 22.0;
  fvec[3] = 27.0;


  answer = new double[n+1];
  answer = solveTridiaognal(n, d, u, l, fvec);

  //Print result
  for(int i=0; i<=n; i++)
  {
    std::cout << "c[" << i << "] = " << answer[i] << "\n";
  }

  delete[] answer;
}



int main(int argc, char* argv[])
{
  double *fvec, *d, *u, *l;

  testExample(d, u, l, fvec);

  delete[] fvec;
  delete[] d;
  delete[] u;
  delete[] l;

  return 0;
}
