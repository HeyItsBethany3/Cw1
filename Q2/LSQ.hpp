#ifndef LSQHEADERDEF
#define LSQHEADERDEF

#include "AbstractFunction.hpp"

#include<string>

/* Class for least squares quadrature using the legendre polynomials on the
interval [-1,1] */

class LSQ {

public:

  // Constructor
  LSQ(const int num, AbstractFunction& aFunction);

  // Destructor
  ~LSQ();

  // Legendre polynomials
  static double phi_1(const double x);
  static double phi_2(const double x);
  static double phi_3(const double x);
  static double phi_4(const double x);
  // static means they do not need an object to be instantiated to be used

// ------ Method 1 for calculating bArray using Gaussian functions -------

  // Find 3-pont Gaussian approximation on interval [-1,1] given function f
  double Gauss3(double (*P_function)(const double x));

  // Find 5-pont Gaussian approximation on interval [-1,1] given function f
  double Gauss5(double (*P_function)(const double x));

  // Finds b values
  void findBGauss3(); // Calculates bArray using 'Gauss3'
  void findBGauss5(); // Calculates bArray using 'Gauss5'

// ----- Method 2 for calculating bArray using different gaussian functions -------

  // Alternative 3-point gaussian method
  double altGauss3(double (*P_function)(const double x));

  // Alternative 5-point gaussian method
  double altGauss5(double (*P_function)(const double x));

  // Finds b values
  void altFindBGauss3(); // Calculates bArray using 'altGauss3'
  void altFindBGauss5(); // Calculates bArray using 'altGauss5'

// -----------------------------------------------------------------------

  // Finds c coefficients (given that b has already been calculated)
  void computeCoefficients();

  void showB(); // Shows bArray vector
  void showC(); // Shows cArray coefficients

  // Evaluate least squares approximation q_n at point x
  double evaluateQ(const double x);

  // Plots q_n values (discretising interval into nodes+1 points)
  void plotQ(const int nodes);

  // Finds an L2 error norm approximation by discretising interval x_k
  double errorNorm(const int nodes); // nodes+1 = number of discretisation points

protected:

  int n;  // Order of LSQ approx polynomial q_n
  double* bArray; // Vector of <f, p> terms where p are the legendre polynomials (length 4)
  double* cArray; // Vector of c coefficients (length n)
  double* cFull; // Complete vector of c coefficients (length 4)
  AbstractFunction* mFunction; // function f(x) to approximate

};

#endif
