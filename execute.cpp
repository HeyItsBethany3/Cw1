#include "Spline.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include<iostream>

int main(int argc, char* argv[]) {

  Function1* f1 = new Function1();
  Spline* s1 = new Spline(10,5,*f1);
  (*s1).Nodes();
  std::cout << "Length " << (*s1).GetLength()<< std::endl;

  delete f1;
  delete s1;

  return 0;
}
