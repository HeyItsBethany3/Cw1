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
  y[1] = y2;
}

int main(int argc, char* argv[]) {
  const double theta0 = M_PI/double(2.0); // initial theta value

  double* yVal;
  yVal = new double[2];
  // TODO: Is my storage allocation right or is this a bad idea?
  euler(theta0,2,8,100, yVal);
  std::cout << "\ny1: " << yVal[0] << " y2: " << yVal[2] << "\n";

  // Testing convergence
  // TODO: What am I looking for? Are y supposed to tend to 0?
  // TODO: Why are y1 and y2 the same?

  double alpha = 2;
  const double T = 7.4162987092/sqrt(alpha); // final time step
  std::cout << "\nT " << T << "\n";
  std::cout << "\nInitial values: y1: " << theta0 << " y2: " << 0;
  int sim = 20;
  int n = 2;
  std::cout.precision(10);
  for(int i=0; i<sim; i++) {
    double* yVal2;
    yVal2 = new double[2];
    euler(theta0,alpha,T, n, yVal2);
    std::cout << "\nh: " << T/double(n)<<" y1: " << yVal2[0] << " y2: " << yVal2[1];
    n= n*2;
    delete[] yVal2;
  }
  std::cout << std::endl;



  delete[] yVal;


  return 0;
}
