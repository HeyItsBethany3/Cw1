#include "LSQ.hpp"
#include "AbstractFunction.hpp"
#include <cmath>
#include <iostream>
#include <cassert>
#include <string>
#include <fstream>

// Constructor
LSQ::LSQ(const int num, AbstractFunction& aFunction) {
  assert(num<=4); // n cannot be greater than 4
  bArray = new double[4]; // Vector of <f, p> terms where p are the legendre polynomials
  cArray = new double[num];  // Vector of c coefficients (length n)
  cFull = new double[4]; // Complete vector of c coefficients (length 4)
  n = num; // n the degree of the polynomial p to find
  mFunction = &aFunction; // function f(x) to approximate
}

// Destructor
LSQ::~LSQ() {
  // Deallocates storage
  delete[] bArray;
  delete[] cArray;
  delete[] cFull;
}

// First legendre polynomial
double LSQ::phi_1(const double x)
{
  return(1.0);
}

// Second legendre polynomial
double LSQ::phi_2(const double x)
{
  return(x);
}

// Third legendre polynomial
double LSQ::phi_3(const double x)
{
  return(0.5*(3.0*pow(x,2)-1));
}

// Fourth legendre polynomial
double LSQ::phi_4(const double x)
{
  return(0.5*(5.0*pow(x,3) - 3.0*x));
}

// Find 3-pont Gaussian approximation on interval [-1,1] given function f
double LSQ::Gauss3(double (*P_function)(const double x))
{
  double sum = (double(8)/double(9))*(*mFunction).evaluateF(0)*(*P_function)(0);
  double xVal = sqrt(double(3)/double(5));
  sum += (double(5)/double(9))*((*mFunction).evaluateF(xVal)*(*P_function)(xVal)+(*mFunction).evaluateF(-xVal)*(*P_function)(-xVal));
  return sum;
}

// Find 5-pont Gaussian approximation on interval [-1,1] given function f
double LSQ::Gauss5(double (*P_function)(const double x)) {
  double sum = (double(128)/double(225))*(*mFunction).evaluateF(0)*(*P_function)(0);
  double xVal = (double(1)/double(3))*sqrt(5.0-(2.0*sqrt(double(10)/double(7))));
  double w1 = double(322+13*double(sqrt(70)))/double(900);
  double xVal2 = (double(1)/double(3))*sqrt(5.0+(2.0*sqrt(double(10)/double(7))));
  double w2 = double(322-13*double(sqrt(70)))/double(900);
  sum += w1*((*mFunction).evaluateF(xVal)*(*P_function)(xVal)+(*mFunction).evaluateF(-xVal)*(*P_function)(-xVal));
  sum += w2*((*mFunction).evaluateF(xVal2)*(*P_function)(xVal2)+(*mFunction).evaluateF(-xVal2)*(*P_function)(-xVal2));
  return sum;
}

// Find B using gauss3 formula
void LSQ::findBGauss3() {
  bArray[0] = LSQ::Gauss3(LSQ::phi_1);
  bArray[1] = LSQ::Gauss3(LSQ::phi_2);
  bArray[2] = LSQ::Gauss3(LSQ::phi_3);
  bArray[3] = LSQ::Gauss3(LSQ::phi_4);
}
/* Note: you can call a class method within another class method but in order
to use a function pointer to another class method, that method must be static
*/

// Find B using gauss5 formula
void LSQ::findBGauss5() {
  bArray[0] = LSQ::Gauss5(LSQ::phi_1);
  bArray[1] = LSQ::Gauss5(LSQ::phi_2);
  bArray[2] = LSQ::Gauss5(LSQ::phi_3);
  bArray[3] = LSQ::Gauss5(LSQ::phi_4);
}


// Alternative 3-point gaussian method
double LSQ::altGauss3(double (*P_function)(const double x)) {
  double x0, x1, x2;

  x0 = sqrt(3) /sqrt(5);
  x1 = 0.;
  x2 = - x0;

  double sum = (8./9.)*(*mFunction).evaluateF(x1)*(*P_function)(x1);
  sum += (5./9.)*(*mFunction).evaluateF(x0)*(*P_function)(x0);
  sum += (5./9.)*(*mFunction).evaluateF(x2)*(*P_function)(x2);
  return sum;

}

// Alternative 5-point gaussian method
double LSQ::altGauss5(double (*P_function)(const double x)) {
  double x0, x1, x2, x3, x4;
  double f0, f1, f2, f3, f4;

  x0 = 0.;
  x1 = (1./3.) * sqrt(5. - 2. * sqrt(10./7.));
  x2 = (-1./3.) * sqrt(5. - 2. * sqrt(10./7.));
  x3 = (1./3.) * sqrt(5. + 2.*sqrt(10./7.));
  x4 = (-1./3.) * sqrt(5. + 2.*sqrt(10./7.));

  f0 = (*mFunction).evaluateF(x0)*(*P_function)(x0);
  f1 = (*mFunction).evaluateF(x1)*(*P_function)(x1);
  f2 = (*mFunction).evaluateF(x2)*(*P_function)(x2);
  f3 = (*mFunction).evaluateF(x3)*(*P_function)(x3);
  f4 = (*mFunction).evaluateF(x4)*(*P_function)(x4);

  double sum = ((128./225.) * f0)+ ((322. + 13. * sqrt(70.))/900.) *(f1+f2);
  sum += ((322. - 13. * sqrt(70.))/900.)*(f3+f4);

  return sum;
}

// Calculates bArray using 'altGauss3'
void LSQ::altFindBGauss3() {
  bArray[0] = LSQ::altGauss3(LSQ::phi_1);
  bArray[1] = LSQ::altGauss3(LSQ::phi_2);
  bArray[2] = LSQ::altGauss3(LSQ::phi_3);
  bArray[3] = LSQ::altGauss3(LSQ::phi_4);
}

// Calculates bArray using 'altGauss5'
void LSQ::altFindBGauss5() {
  bArray[0] = LSQ::altGauss5(LSQ::phi_1);
  bArray[1] = LSQ::altGauss5(LSQ::phi_2);
  bArray[2] = LSQ::altGauss5(LSQ::phi_3);
  bArray[3] = LSQ::altGauss5(LSQ::phi_4);
}


// Finds coefficients of the LSQ approx polynomial q_n
void LSQ::computeCoefficients() {
  double *Mu_inv_diag; //Vector holding the diagonal values of Mu inverse matrix
  Mu_inv_diag = new double[4];

  // Values were found analytically
  Mu_inv_diag[0] = 0.5;
  Mu_inv_diag[1] = 1.5;
  Mu_inv_diag[2] = 2.5;
  Mu_inv_diag[3] =double(7)/double(2);

  for(int i=0; i<4; i++) {
    cFull[i] = Mu_inv_diag[i] * bArray[i]; // Finds all c values
  }
  for(int i=0; i<n; i++)
  {
    cArray[i] = cFull[i]; // Only necessary c values are used
  }

  delete[] Mu_inv_diag;
}

// Displays bArray vector
void LSQ::showB() {
  for(int i=0; i<n; i++)
  {
    std::cout << "b[" << i+1 << "] = " << bArray[i] << "\n";
  }
}

// Displays coefficients cArray
void LSQ::showC() {
  for(int i=0; i<n; i++)
  {
    std::cout << "c[" << i+1 << "] = " << cArray[i] << "\n";
  }
}

// Evaluates q_n at point x on interval [-1,1]
double LSQ::evaluateQ(const double x) {

  double* product = new double[4];
  // All products are calculated
  product[0] = phi_1(x)*cFull[0];
  product[1] = phi_2(x)*cFull[1];
  product[2] = phi_3(x)*cFull[2];
  product[3] = phi_4(x)*cFull[3];

  // Only necessary products are used depending on n
  double sum = 0;
  for(int i=0; i<n; i++) {
    sum += product[i];
  }

  delete[] product;
  return sum;
}

/* Outputs 'nodes', q and f to a file, which is then used to plot q_n against
f(x) using a matlab script 'plotQ.m'.
*/
void LSQ::plotQ(const int nodes) {

  // Writes to file
  std::ofstream file;
  file.open("q.csv");
  assert(file.is_open()); // Checks file opened correctly

  double x, q, exact, error;
  file << nodes << ","<<std::endl; // Outputs number of discretisation points
  double* fvec = new double[nodes+1];

  // Creates nodes on [-1,1] (to evaluate p_n and f at)
  for (int i=0; i<=nodes; i++) {
    x =-1.0 +i*(2.0/double(nodes));
    q = LSQ::evaluateQ(x);
    fvec[i] = (*mFunction).evaluateF(x);
    file << q << ","; // Outputs q values
  }
  file << std::endl;

  for (int i=0; i<=nodes; i++) {
    file << fvec[i] << ","; // Outputs f values
  }
  delete[] fvec;

  // Closes file and moves file to correct directory for matlab file to be easily run
  file.close();
  std::system("rm ../../../MATLAB/q.csv");
  std::system("cp q.csv ../../../MATLAB");
}

// Approximates the L2 norm of q_n, by discretising the interval
double LSQ::errorNorm(const int nodes) {

  double x, q, exact, error;
  double norm = 0;

  // Discretises the interval [-1,1] into nodes+1 points
  for (int i=0; i<=nodes; i++) {
    x =-1.0 +i*(2.0/double(nodes));
    q = LSQ::evaluateQ(x);
    exact = (*mFunction).evaluateF(x);
    error = fabs(q-exact); // Finds the error at each point
    norm += pow(error, 2);

  }
  return sqrt(norm); // Outputs L2 norm 
}
