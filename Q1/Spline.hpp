#ifndef SPLINEHEADERDEF
#define SPLINEHEADERDEF

#include "AbstractFunction.hpp"

class Spline {

  public:
    // Constructor
    Spline(const double len, const int n, AbstractFunction& aFunction);

    // Destructor
    ~Spline();

    // Construct nodes
    void Nodes();

    // Find system of equations
    void FindSystem();

    // Shows system of equations
    void showSystem();

    // Solves system of equations (finds coefficients of spline)
    void solveTridiaognal();

    // shows coefficients of spline
    void showCoeff();

  protected:

    double mLen; // Interval is (0, mLen)
    int mN; // n+1 interpolating nodes
    double *mNodes; // x nodes
    double mH;

    double *mFvec;
    double *mDiag;
    double *mUpper;
    double *mLower;

    double *mCoeff;

    AbstractFunction* mFunction;


};


#endif
