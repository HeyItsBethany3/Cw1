#include "LSQ.hpp"
#include "Function1.hpp"
#include "Function2.hpp"
#include "AbstractFunction.hpp"
#include <iostream>

int main(int argc, char* argv[]) {


  // Part c and d

  Function1* f1 = new Function1(); // Initialises function object
  LSQ* lsq =  new LSQ(4, *f1); // Initialises least squares problem

  // Using 3-point gauss quadrature
  std::cout <<"\n3-point Gauss quadrature:\n";
  (*lsq).findBGauss3();
  (*lsq).computeCoefficients();
  //(*lsq).showB();
  (*lsq).showC();


  std::cout << "\nq(0.2):\t" << (*lsq).evaluateQ(0.2)<< std::endl;
  std::cout << "f(0.2):\t" << (*f1).evaluateF(0.2)<< std::endl;
  (*lsq).plotQ(100);


  // Using 5-point gauss quadrature
  std::cout <<"\n5-point Gauss quadrature:\n";
  LSQ* lsq2 =  new LSQ(4, *f1);
  (*lsq2).findBGauss5();
  (*lsq2).computeCoefficients();
  (*lsq2).showC();

  std::cout << "\nq(0.2):\t" << (*lsq2).evaluateQ(0.2)<< std::endl;
  std::cout << "f(0.2):\t" << (*f1).evaluateF(0.2)<< std::endl;
  std::cout << "Error norm: " << (*lsq2).errorNorm(100)<< std::endl;


  // Same thing repeated using simple function 1+x+x^2
  //Function* f2 = new Function2();

  // Part e


  // Calculate error norms for n=1,2,3,4
  int sim = 100;
  for(int i=1; i<=4; i++) {
    Function2* f2 = new Function2(); // Initialises function object
    LSQ* lsq3 =  new LSQ(i, *f2); //
    (*lsq3).findBGauss5();
    (*lsq3).computeCoefficients();
    std::cout << "\nError norm for q_" << i << " : " << (*lsq3).errorNorm(sim);

    delete lsq3;
    delete f2;
  }
  std::cout << std::endl;


  delete lsq;
  delete f1;
  delete lsq2;

  return 0;
}
