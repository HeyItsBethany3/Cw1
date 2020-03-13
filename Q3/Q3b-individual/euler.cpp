#include <iostream>
#include <cmath>

/* Implements euler's method
theta0: initial theta value, T is the time step of interest, alpha is a parameter
of the problem and n is the number of time intervals, y is an empty vector of length 2
Finds vector of (y1,y2) approximate values at time T */
void euler(const double theta0, const double alpha,const double T, const int n, double* y) {

  const double h = T/double(n); // step-size
  double y1old = theta0;
  double y2old = 0;
  double y1, y2;
  for (int i=1; i<=n; i++) {
    y1 = y1old + (h*y2old);
    y2 = y2old + (-h*pow(alpha,2)*sin(y1old));
    y1old = y1;
    y2old = y2;
  }

  y[0] = y1;
  y[1] = y[0];
}

int main(int argc, char* argv[]) {
  const double theta0 = M_PI/double(2.0); // initial theta value

  double* yVal;
  yVal = new double[2];
  // TODO: Is my storage allocation right or is this a bad idea?
  euler(theta0,2,8, 100, yVal);
  std::cout << "\ny1: " << yVal[0] << " y2: " << yVal[2] << "\n";

  // Testing convergence
  double alpha = 2;
  const double T = 7.4162987092/sqrt(alpha); // final time step
  int sim = 10;
  int n = 10;
  std::cout.precision(10);
  for(int i=0; i<sim; i++) {
    euler(theta0,alpha,T, 100, yVal);
    std::cout << "\nh: " << T/double(n)<<" y1: " << yVal[0] << " y2: " << yVal[2] << "\n";
    n= n*2;
  }



  delete[] yVal;


  return 0;
}
