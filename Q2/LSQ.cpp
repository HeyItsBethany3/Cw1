#include "LSQ.hpp"
#include "AbstractFunction.hpp"
#include <cmath>
#include <iostream>
#include <cassert>


LSQ::LSQ(const int num, AbstractFunction& aFunction) {
  assert(num<=4); // n cannot be greater than 4
  bArray = new double[4];
  cArray = new double[num];
  n = num;
  mFunction = &aFunction;
}

LSQ::~LSQ() {

  delete[] bArray;
  delete[] cArray;
}

double LSQ::phi_1(const double x)
{
  return(1.0);
}

double LSQ::phi_2(const double x)
{
  return(x);
}

double LSQ::phi_3(const double x)
{
  return(0.5*(3.0*pow(x,2.0)-1));
}

double LSQ::phi_4(const double x)
{
  return(0.5*(5.0*pow(x,3.0) - 3.0*x));
}



double LSQ::Gauss3(double (*P_function)(const double x))
{
  double sum = (double(8)/double(9))*(*mFunction).evaluateF(0)*(*P_function)(0);
  double xVal = sqrt(double(3)/double(5));
  sum += (double(5)/double(9))*((*mFunction).evaluateF(xVal)*(*P_function)(xVal)+(*mFunction).evaluateF(-xVal)*(*P_function)(-xVal));
  return sum;
}


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

// Finds coefficients of the LSQ approx polynomial q_n
void LSQ::computeCoefficients() {
  double *Mu_inv_diag; //Vector holding the diagonal values of Mu inverse matrix

  Mu_inv_diag = new double[4];

  // Values were found analytically
  Mu_inv_diag[0] = 0.5;
  Mu_inv_diag[1] = 1.5;
  Mu_inv_diag[2] = 2.5;
  Mu_inv_diag[3] = 9.0;

  for(int i=0; i<n; i++)
  {
    cArray[i] = Mu_inv_diag[i] * bArray[i];
  }

  delete[] Mu_inv_diag;
}

void LSQ::showB() {
  for(int i=0; i<n; i++)
  {
    std::cout << "b[" << i+1 << "] = " << bArray[i] << "\n";
  }
}

void LSQ::showC() {
  for(int i=0; i<n; i++)
  {
    std::cout << "c[" << i+1 << "] = " << cArray[i] << "\n";
  }
}
