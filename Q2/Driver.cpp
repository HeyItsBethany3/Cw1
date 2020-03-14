#include "LSQ.hpp"
#include "Function1.hpp"
#include "Function2.hpp"
#include "AbstractFunction.hpp"
#include <iostream>
#include <string>

// quadrature: gauss method to use, n is order to approximate
void showCoefficients(AbstractFunction* funct, const int n,  std::string quadrature) {
    LSQ* lsq =  new LSQ(n, *funct); // Initialises least squares problem
    std::cout << "\nLeast squares approximation q_" << n << " using quadrature ";
    std::cout << quadrature  << "\n";
    if (quadrature=="Gauss3") {
        (*lsq).findBGauss3();
    }
    else if (quadrature == "Gauss5") {
        (*lsq).findBGauss5();
    }
    else if (quadrature == "altGauss3") {
        (*lsq).altFindBGauss3();
    }
    else if (quadrature == "altGauss5") {
        (*lsq).altFindBGauss5();
    }

    (*lsq).computeCoefficients();
    (*lsq).showC();

    delete lsq;
}

// Calculates error norms for n=1,2,3,4, sim is the number of discretisation points
void showError(AbstractFunction* funct, const int sim, std::string quadrature) {
  for(int i=1; i<=4; i++) {
    LSQ* lsq3 =  new LSQ(i, *funct); //
    if (quadrature=="Gauss5") {
        (*lsq3).findBGauss5();
    }
    else if (quadrature == "altGauss5") {
        (*lsq3).altFindBGauss5();
    }

    (*lsq3).computeCoefficients();
    std::cout << "\nError norm for q_" << i << " : " << (*lsq3).errorNorm(sim);

    delete lsq3;

  }
  std::cout << std::endl;

}

// Function prototypes
void showCoefficients(AbstractFunction* funct, const int n,  std::string quadrature);
void showError(AbstractFunction* funct, const int sim, std::string quadrature);

int main(int argc, char* argv[]) {

  // Part c and d: for function in question
  Function1* f1 = new Function1(); // Initialises function object
  showCoefficients(f1, 4, "Gauss3");   // Using 3-point gauss quadrature
  showCoefficients(f1, 4, "Gauss5");   // Using 5-point gauss quadrature
  showCoefficients(f1, 4, "altGauss3");   // 3-point alternative gauss method
  showCoefficients(f1, 4, "altGauss5");   // 5-point alternative gauss method

  // Part e
  // Calculate error norms for n=1,2,3,4
  showError(f1, 100, "Gauss5");
  showError(f1, 100, "altGauss5");
  /* We get high errors but this is because the function f(x) is hard to
  approximate using low order polynomials */

  // Plotting q_4 and f(x) (using matlab script)
  LSQ* lsq2 =  new LSQ(4, *f1);
  (*lsq2).findBGauss5();
  (*lsq2).computeCoefficients();
  (*lsq2).plotQ(100); // plots Q using 100 discretisation points

  // Repeating part e for the simple function f(x)=1+x+x^2 to prove that the code works
  Function2* f2 = new Function2();
  showError(f2, 100, "Gauss5");

  delete f1;
  delete f2;

  return 0;
}
