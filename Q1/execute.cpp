#include "Spline.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include<iostream>




int main(int argc, char* argv[]) {

// Question 1c and 1d
  int n = 8;
  double len = 2.0;

  Function1* f1 = new Function1();
  Spline* s1 = new Spline(len,n,*f1);
  (*s1).Nodes();
  (*s1).FindSystem();
  (*s1).showSystem();

  (*s1).solveTridiaognal();
  (*s1).showCoeff();


  delete f1;
  delete s1;

  return 0;
}
