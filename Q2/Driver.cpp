#include "LSQ.hpp"
#include <iostream>

int main(int argc, char* argv[]) {


  LSQ* obj1 =  new LSQ(4);
  std::cout <<"\nTest " << (*obj1).phi_2(0.3) << std::endl;



  delete obj1;

  return 0;
}
