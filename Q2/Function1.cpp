#include "Function1.hpp"
#include<cmath>

// Constructor
Function1::Function1() {
}

// Specifies f(x)
double Function1::evaluateF(double x) {
  return (exp(pow(x,2)+x) * sin((double(5)/double(4))*M_PI*x));

}


// Destructor
Function1::~Function1() {

}
