#include "Spline.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include<iostream>


int main(int argc, char* argv[]) {

// Question 1c and 1d
  int n = 8;
  double len = 2.0;

  Function1* f1 = new Function1(); // Initialises function object
  Spline* s1 = new Spline(len,n,*f1); // Initialises spline object

  (*s1).Nodes(); // Creates nodes
  (*s1).FindSystem(); // Finds system of equations
  (*s1).showSystem();

  (*s1).solveTridiaognal(); // Solves system of equations
  (*s1).showCoeff();


  delete f1;
  delete s1;

  return 0;
}
