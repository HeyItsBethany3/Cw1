#include "Function2.hpp"
#include<cmath>

// Constructor
Function2::Function2() {
}

// Specifies f(x)
double Function2::evaluateF(double x) {
  return 1+x+pow(x,2);
}


// Destructor
Function2::~Function2() {

}
