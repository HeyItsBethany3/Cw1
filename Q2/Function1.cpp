#include "Function1.hpp"
#include<cmath>

// Constructor
Function1::Function1() {
}

// Specifies f(x)
double Function1::evaluateF(double x) {
  return exp(pow(x,2.0)+x) * sin(1.25*M_PI*x);
}


// Destructor
Function1::~Function1() {

}
