#ifndef SPLINEHEADERDEF
#define SPLINEHEADERDEF

#include "AbstractFunction.hpp"

// Class for constructing a spline to aprroximate a function over [0,len]
class Spline {

  public:
    // Constructor: spline approximates 'aFunction' over interval [0,len] with order n
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

    // Alternative method for solving the system
    void solveMethod2(); // - Ioanna

    // Shows coefficients of spline
    void showCoeff();

    // Finds B (cubic spline)
    double evaluateB(const double xstar);

    // Evaluate spline at x
    double evaluateSpline(const double x);

    // Find absolute error of spline approximation to f at point x
    double error(const double x);


  protected:

    double mLen; // Interval is (0, mLen)
    int mN; // mN+1 interpolating nodes
    double *mNodes; // x_0, x_1, ... , x_n
    double *mFullNodes; // Complete set of nodes x_(-1), x_0, x_1, ... , x_n, x_(n+1)
    double mH; // Distance between nodes

    double *mFvec; // f(x) at all nodes x
    double *mDiag; // Diagonal
    double *mUpper; // Upper diagonal
    double *mLower; // Lower diagonal

    double *mCoeff; // Spline coefficients
    double *mFullC; // Complete set of spline coefficients

    AbstractFunction* mFunction; // Function pointer to f

};

#endif
