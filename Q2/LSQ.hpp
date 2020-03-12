#ifndef LSQHEADERDEF
#define LSQHEADERDEF

#include "AbstractFunction.hpp"

/* Class for least squares quadrature using the legendre polynomials on the
interval [-1,1] */

class LSQ {

public:

  // Constructor
  LSQ(const int num, AbstractFunction& aFunction);

  // Destructor
  ~LSQ();

  // Legendre polynomials
  static double phi_1(const double x); // static means they do not need an object to be instantiated to be used
  static double phi_2(const double x);
  static double phi_3(const double x);
  static double phi_4(const double x);

  // Find b values
  void findBGauss3();
  void findBGauss5();

  // Find 3-pont Gaussian approximation on interval [-1,1] given function f
  double Gauss3(double (*P_function)(const double x));

  // Find 5-pont Gaussian approximation on interval [-1,1] given function f
  double Gauss5(double (*P_function)(const double x));

  void computeCoefficients();


protected:

  int n;  // Order of LSQ approx polynomial q_n

  double* bArray; // Vector of <f, p> terms
  double* cArray; // Vector of c coefficients


  AbstractFunction* mFunction;


};

#endif
