#include <iostream>
#include <cmath>

/* Implements euler's method
theta0 is initial y1 value, T is the time step of interest, alpha is a parameter
of the problem and n is the number of time intervals, y is an empty vector of length 2
Finds vector of (y1,y2) approximate values at time T */
void euler(const double theta0, const double alpha,const double T, const int n, double* y) {

  const double h = T/double(n); // step-size
  double y1old = theta0; // Initial y1 value
  double y2old = 0; // Initial y2 value
  double y1, y2;
  for (int i=1; i<=n; i++) {
    // Finds next y1 and y2 value using euler's method
    y1 = y1old + (h*y2old);
    y2 = y2old + (-h*pow(alpha,2)*sin(y1old));
    y1old = y1;
    y2old = y2;
  }
  // Retrieves y1 and y2 approximations at time T (after n steps)
  y[0] = y1;
  y[1] = y2;
}

// Function prototypes
void euler(const double theta0, const double alpha,const double T, const int n, double* y);

int main(int argc, char* argv[]) {

  // Simple test
  const double theta0 = M_PI/double(2.0); // initial theta value
  double* yVal;
  yVal = new double[2];
  euler(theta0,2,8,100,yVal);
  std::cout << "\ny1: " << yVal[0] << " y2: " << yVal[2] << "\n";

  // Testing convergence
  double alpha = 2;
  const double T = 7.4162987092/(alpha); // final time step
  std::cout << "\nInitial values: y1: " << theta0 << " y2: " << 0;
  int sim = 20; // number of times we run test for
  int n = 2; // initial n value
  std::cout.precision(10);

  // With these parameters we expect the results to be y1 = pi/2 and y2=0
  for(int i=0; i<sim; i++) {
    double* yVal2;
    yVal2 = new double[2];
    euler(theta0,alpha,T, n, yVal2);
   
    // Output approximations found
    std::cout << "\nh: " << T/double(n)<<"\ty1: " << yVal2[0] << "\ty2: " << yVal2[1];
    // Double n each time (halve h)
    n=n*2;

    delete[] yVal2;
  }
    std::cout << std::endl;

  // As h decreases (it halves each time), the y1 and y2 values converge
  // They also converge to the expected values!

  delete[] yVal;

  return 0;
}
