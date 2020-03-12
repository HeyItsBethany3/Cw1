#include "AbstractFunction.hpp"

class Function1: public AbstractFunction {

  public:
    // Constructor
    Function1();

    // Destructor
    ~Function1();

    // f(x)
    double evaluateF(double x);


};
