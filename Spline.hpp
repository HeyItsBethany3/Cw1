#ifndef SPLINEHEADERDEF
#define SPLINEHEADERDEF

#include "AbstractFunction.hpp"

class Spline {

  public:
    // Constructor
    Spline(const double len, const int n, AbstractFunction& aFunction);

    // Destructor
    ~Spline();

    // Get length
    double GetLength();

    // Construct nodes
    void Nodes();

    // Find system of equations
    void FindSystem();

  protected:

    double mLen; // Interval is (0, mLen)
    int mN; // n+1 interpolating nodes
    double *mNodes; // x nodes
    double mH;

    double *mFvec;
    double *mDiag;
    double *mUpper;
    double *mLower;

    AbstractFunction* mFunction;


};


#endif
