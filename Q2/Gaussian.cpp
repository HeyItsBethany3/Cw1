#include<cmath>
#include<iostream>
/* Implement 3-point and 5-point Gaussian quadrature to approximate integrals
on [-1,1].
*/

// Monomial testing function
double f(const double x, const int n) {
  return pow(x,n);
}

// 3-point Gaussian quadrature rule
double Gauss3(double (*function)(const double x, const int n), const int n) {
  double sum = (double(8)/double(9))*(*function)(0, n);
  double xVal = sqrt(double(3)/double(5));
  sum += (double(5)/double(9))*((*function)(xVal,n)+(*function)(-xVal,n));
  return sum;
}

int main(int argc, char* argv[]) {

  double val = Gauss3(f,2);
  std::cout << val;





  return 0;
}
