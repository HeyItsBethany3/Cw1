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

  // ----- Alternative methods ------
  // Alternative 3-point gaussian method
  double altGauss3(double (*P_function)(const double x));
  void altFindBGauss3();

  // ---------------------------------


  // Finds c coefficients (given that b has already been calculated)
  void computeCoefficients();

  void showB(); // shows b vector

  void showC(); // shows c coefficients

  // evaluate least squares approximation q_n at point x
  double evaluateQ(const double x);

  // Finds an L2 error norm approximation by discretising interval x_k
  double errorNorm(const int nodes); // nodes+1 = number of discretisation points

  // plots q values (discretising interval )
  void plotQ(const int nodes);


protected:

  int n;  // Order of LSQ approx polynomial q_n

  double* bArray; // Vector of <f, p> terms (length 4)
  double* cArray; // Vector of c coefficients (length n)
  double* cFull; // Complete vector of c coefficients (length 4)
  AbstractFunction* mFunction;


};

#endif
