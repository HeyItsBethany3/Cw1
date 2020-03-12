#include "LSQ.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include <iostream>



int main(int argc, char* argv[]) {


  Function1* f1 = new Function1(); // Initialises function object
  LSQ* obj1 =  new LSQ(4, *f1); // Initialises least squares problem
  std::cout <<"\nTest " << (*obj1).phi_2(0.3) << std::endl;

  std::cout << (*f1).evaluateF(0.3) << std::endl;
  (*obj1).findBGauss3();




  delete obj1;
  delete f1;

  return 0;
}
