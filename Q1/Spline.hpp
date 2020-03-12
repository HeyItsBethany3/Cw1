#ifndef SPLINEHEADERDEF
#define SPLINEHEADERDEF

#include "AbstractFunction.hpp"

class Spline {

  public:
    // Constructor
    Spline(const double len, const int n, AbstractFunction& aFunction);

    // Destructor
    ~Spline();

    // Constructs nodes
    void Nodes();

    // Finds system of equations
    void FindSystem();

    // Shows system of equations
    void showSystem();

    // Solves system of equations (finds coefficients of spline)
    void solveTridiaognal(); // - Britta

    void solveMethod2(); // - Ioanna

    double evaluateB(const double xstar); // Finds B

    // Shows coefficients of spline
    void showCoeff();

    // Evaluate spline at x
    double evaluateSpline(const double x);


    // Find absolute error of spline approximation to f at point x
    double error(const double x);




  protected:

    double mLen; // Interval is (0, mLen)
    int mN; // mN+1 interpolating nodes
    double *mNodes;
    double *mFullNodes; //complete set of nodes
    double mH; // Distance between nodes

    double *mFvec; // f(x) at all nodes x
    double *mDiag; // Diagonal
    double *mUpper; // Upper diagonal
    double *mLower; // Lower diagonal

    double *mCoeff; // Spline coefficients
    double *mFullC; // All spline coefficients

    AbstractFunction* mFunction; // Function pointer to f

};

#endif
