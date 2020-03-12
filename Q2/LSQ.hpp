#ifndef LSQHEADERDEF
#define LSQHEADERDEF

/* Class for least squares quadrature using the legendre polynomials on the
interval [-1,1] */

class LSQ {

public:

  // Constructor
  LSQ(const int num);

  // Destructor
  ~LSQ();

  // Legendre polynomials
  double phi_1(const double x);
  double phi_2(const double x);
  double phi_3(const double x);
  double phi_4(const double x);

protected:

  int n;  // Order of LSQ approx polynomial q_n

  double* bArray; // Vector of <f, p> terms
  double* cArray; // Vector of c coefficients




};

#endif
