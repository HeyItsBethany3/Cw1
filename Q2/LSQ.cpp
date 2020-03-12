#include "LSQ.hpp"
#include "AbstractFunction.hpp"
#include <cmath>
#include <iostream>


LSQ::LSQ(const int num, AbstractFunction& aFunction) {
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


// Find B using gauss5 formula
void LSQ::findBGauss5() {
  bArray[0] = LSQ::Gauss5(LSQ::phi_1);
  bArray[1] = LSQ::Gauss5(LSQ::phi_2);
  bArray[2] = LSQ::Gauss5(LSQ::phi_3);
  bArray[3] = LSQ::Gauss5(LSQ::phi_4);
}

void LSQ::computeCoefficients() {

}
