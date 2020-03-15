#include<cmath>
#include<iostream>
/* Implement 3-point and 5-point Gaussian quadrature to approximate integrals
on [-1,1].
Both Gaussian methods are designed for a monomial function, so would have to
be modified for other functions f(x)
*/

// Monomial testing function
double f(const double x, const int n) {
  return pow(x,n);
}

// 3-point Gaussian quadrature rule
double Gauss3(double (*function)(const double x, const int n), const int n) {
  /* Takes function pointer as input and n which is the degree of the power in f
   Outputs integral approximation */
  double sum = (double(8)/double(9))*(*function)(0, n);
  double xVal = sqrt(double(3)/double(5));
  sum += (double(5)/double(9))*((*function)(xVal,n)+(*function)(-xVal,n));
  return sum;
}

// 5-point Gaussian quadrature rule
double Gauss5(double (*function)(const double x, const int n), const int n) {
  /* Takes function pointer as input and n which is the degree of the power in f
   Outputs integral approximation */

  double sum = (double(128)/double(225))*(*function)(0, n);
  double xVal = (double(1)/double(3))*sqrt(5.0-(2.0*sqrt(double(10)/double(7))));
  double w1 = double(322+13*double(sqrt(70)))/double(900);
  double xVal2 = (double(1)/double(3))*sqrt(5.0+(2.0*sqrt(double(10)/double(7))));
  double w2 = double(322-13*double(sqrt(70)))/double(900);
  sum += w1*((*function)(xVal,n)+(*function)(-xVal,n));
  sum += w2*((*function)(xVal2,n)+(*function)(-xVal2,n));

  return sum;

}

// Function prototypes
double f(const double x, const int n);
double Gauss3(double (*function)(const double x, const int n), const int n);
double Gauss5(double (*function)(const double x, const int n), const int n);

int main(int argc, char* argv[]) {

  // Tests each Gaussian function for different monomial functions
  std::cout.precision(10);
  for (int k=0; k<=10; k++){
    std::cout << "\nk: " << k;
    std::cout << "\t3-point rule: " << Gauss3(f,k);
    std::cout << "\t5-point rule: " << Gauss5(f,k);
  }
  std::cout << std::endl;
  return 0;
}
