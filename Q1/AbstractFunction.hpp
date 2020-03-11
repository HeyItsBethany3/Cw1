#ifndef ABSTRACTFUNCTION
#define ABSTRACTFUNCTION

class AbstractFunction {
  public:
    // f(x)
    virtual double evaluateF(double x) = 0;

    // f'(x)
    virtual double derivative(double x) = 0;

    // Destructor
    virtual ~AbstractFunction() = 0;

    

};

#endif
