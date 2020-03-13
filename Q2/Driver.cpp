#include "LSQ.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include <iostream>



int main(int argc, char* argv[]) {


  Function1* f1 = new Function1(); // Initialises function object
  LSQ* lsq =  new LSQ(3, *f1); // Initialises least squares problem

  // Using 3-point gauss quadrature
  std::cout <<"\n3-point Gauss quadrature:\n";
  (*lsq).findBGauss3();
  (*lsq).computeCoefficients();
  //(*lsq).showB();
  (*lsq).showC();

  // Using 5-point gauss quadrature
  std::cout <<"\n5-point Gauss quadrature:\n";
  (*lsq).findBGauss5();
  (*lsq).computeCoefficients();
  (*lsq).showC();



  delete lsq;
  delete f1;

  return 0;
}
