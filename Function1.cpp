#include "Function1.hpp"
#include<cmath>

Function1::Function1() {
}

double Function1::evaluateF(double x) {
  return (exp(x)* cos(1.25*M_PI*x));
}

Function1::~Function1() {

}

double Function1::derivative(double x) {
  return (exp(x)*(cos(1.25*M_PI*x)-(1.25*M_PI*sin(1.25*M_PI*x))));
}
