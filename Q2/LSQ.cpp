#include "LSQ.hpp"
#include <cmath>


LSQ::LSQ(const int num) {
  bArray = new double[num];
  cArray = new double[num];
  n = num;

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
