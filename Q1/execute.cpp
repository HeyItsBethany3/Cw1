#include "Spline.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include<iostream>


int main(int argc, char* argv[]) {

  // Question 1c and 1d
  int n = 8;
  double len = 2.0; // So h=1/4

  Function1* f1 = new Function1(); // Initialises function object
  Spline* s1 = new Spline(len,n,*f1); // Initialises spline object

  (*s1).Nodes(); // Creates nodes
  (*s1).FindSystem(); // Finds system of equations
  //(*s1).showSystem();

  (*s1).solveTridiaognal(); // Solves system of equations
  (*s1).showCoeff();

  // Question 1e
  std::cout << "\nSpline at x=0.3 is " << (*s1).evaluateSpline(0.3);
  std::cout << "\nf(0.3) is " << (*f1).evaluateF(0.3) << std::endl;

  std::cout << "\nSpline at x=0.9 is " << (*s1).evaluateSpline(0.9);
  std::cout << "\nf(0.9) is " << (*f1).evaluateF(0.9) << std::endl;

  delete f1;
  delete s1;

  return 0;
}
