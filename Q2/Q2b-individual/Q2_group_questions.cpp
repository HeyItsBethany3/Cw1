#include <iostream>
#include <cmath>

#include "Q2_Option2_LSQ_coeffs.cpp" //TODO can't include .cpp

//Note: Gauss functions were slightly adjusted to take 2 functions as input

// 3-point Gaussian quadrature rule
double Gauss3(double (*F_function)(const double x), double (*P_function)(const double x))
{
  double sum = (double(8)/double(9))*(*F_function)(0)*(*P_function)(0);
  double xVal = sqrt(double(3)/double(5));
  sum += (double(5)/double(9))*((*F_function)(xVal)*(*P_function)(xVal)+(*F_function)(-xVal)*(*P_function)(-xVal));
  return sum;
}

// 5-point Gaussian quadrature rule
double Gauss5(double (*F_function)(const double x), double (*P_function)(const double x))
{
  double sum = (double(128)/double(225))*(*F_function)(0)*(*P_function)(0);
  double xVal = (double(1)/double(3))*sqrt(5.0-(2.0*sqrt(double(10)/double(7))));
  double w1 = double(322+13*double(sqrt(70)))/double(900);
  double xVal2 = (double(1)/double(3))*sqrt(5.0+(2.0*sqrt(double(10)/double(7))));
  double w2 = double(322-13*double(sqrt(70)))/double(900);
  sum += w1*((*F_function)(xVal)*(*P_function)(xVal)+(*F_function)(-xVal)*(*P_function)(-xVal));
  sum += w2*((*F_function)(xVal2)*(*P_function)(xVal2)+(*F_function)(-xVal2)*(*P_function)(-xVal2));
  return sum;
}

// TODO: Should we not just have another function which multiplies two functionds together?


double evaluateF(const double x)
{
  return(exp(pow(x,2.0)+x) * sin(1.25*M_PI*x)); //TODO double check written right
}

double* calculateB(double (*F_function)(const double x), double (*Gauss_function)(double (*F_function)(const double x), double (*P_function)(const double x)))
{
  double *b;
  b = new double[4];

  b[0] = (*Gauss_function)((*F_function), phi_1);
  b[1] = (*Gauss_function)((*F_function), phi_2);
  b[2] = (*Gauss_function)((*F_function), phi_3);
  b[3] = (*Gauss_function)((*F_function), phi_4);

  return b;
}

double* computeCoefficients(const int n, double (*F_function)(const double x), double (*Gauss_function)(double (*F_function)(const double x), double (*P_function)(const double x)))
/*
Parameters: n = order of LSQ approx polynomial q_n
*/
{
  double *b, *coeffs;
  b = new double[n];
  b = calculateB((*F_function), (*Gauss_function));

  coeffs = new double[n];
  coeffs = solveLSQCoefficients(n, b);

  //Only outputs the n coefficients for q_n
  for(int i=0; i<n; i++)
  {
    std::cout << "c[" << i+1 << "] = " << coeffs[i] << "\n";
  }

  delete[] b;

  return coeffs;
}


int main(int argc, char* argv[])
{
  std::cout << "3-point Gauss Quadrature:" <<"\n";
  computeCoefficients(4, evaluateF, Gauss3);

  std::cout << "\n\n5-point Gauss Quadrature:" <<"\n";
  computeCoefficients(4, evaluateF, Gauss5);

  return 0;
}
