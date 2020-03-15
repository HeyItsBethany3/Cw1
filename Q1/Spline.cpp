#include "Spline.hpp"
#include "AbstractFunction.hpp"
#include<iostream>
#include<cassert>
#include<cmath>


// Constructor
Spline::Spline(const double len, const int n, AbstractFunction& aFunction) {

  // Checks interval and n is reasonable
  assert(len>0);
  assert(n>0);

  mLen = len; // spline interval is [0,len]
  mN = n; // n+1 interpolating points
  mNodes = new double[n+1]; // stores interpolation points x_0, x_1, ... , x_n
  mFullNodes = new double[n+3]; // stores end nodes x_(-1), x_(n+1) too
  mH = len/double(n); // distance between nodes

  mFvec = new double[n+1]; // f values at each node
  mDiag = new double[n+1]; // diagonal vector
  mUpper = new double[n]; // upper diagonal vector
  mLower = new double[n]; // lower diagonal vector
  mCoeff = new double[n+1]; // coefficients of spline c_0, c_1, ... , c_n
  mFullC = new double[n+3]; // includes end coefficients c_(-1), c_(n+1)

  mFunction = &aFunction; // f function for nodes to be evaluated at

}

// Destructor
Spline::~Spline() {
  // Deallocates storage
  delete mNodes;
  delete mFvec;
  delete mDiag;
  delete mUpper;
  delete mLower;
  delete mCoeff;
  delete mFullC;
  delete mFullNodes;
}

// Constructs interpolating nodes
void Spline::Nodes() {
  for (int i=0; i<mN+1; i++) {
    mNodes[i] = i*mH;
  }

  for (int i=1; i<mN+2; i++) {
    mFullNodes[i] = (i-1)*mH;
  }
  mFullNodes[0] = -mH;
  mFullNodes[mN+2] = mLen +mH;
}

// Finds system of equations for spline
void Spline::FindSystem() {
  // We are creating a system of the form Ac=b where A is a tridiagonal matrix,
  // b is mFvec and c is mCoeff

  // Find the RHS vector for the system of equations, which we call mFvec
  for (int i=1; i<mN; i++) {
    mFvec[i]=(*mFunction).evaluateF(mNodes[i]);
  }
  mFvec[0]=(*mFunction).evaluateF(mNodes[0])+(((double(1)/double(3))*mH*(*mFunction).derivative(mNodes[0])));
  mFvec[mN]=(*mFunction).evaluateF(mNodes[mN])-(((double(1)/double(3))*mH*(*mFunction).derivative(mNodes[mN])));

  // Find diagonal elements of A
  for (int i=0; i<=mN; i++) {
    mDiag[i] = 4;
  }

  // Construct upper diagonal elements of A
  mUpper[0]=2;
  for (int i=1; i<mN; i++) {
    mUpper[i]=1;
  }

  // Construct lower diagonal elements of A
  mLower[mN-1]=2;
  for (int i=0; i<mN-1; i++) {
    mLower[i]=1;
  }
}

// Displays system to solve
void Spline::showSystem() {
  std::cout << "\nF(x): ";
  for (int i=0; i<=mN; i++) {
    std::cout << mFvec[i] << " ";
  }

  std::cout << "\nd: ";
  for (int i=0; i<=mN; i++) {
    std::cout << mDiag[i] << " ";
  }

  std::cout << "\nu: ";
  for (int i=0; i<mN; i++) {
    std::cout << mUpper[i] << " ";
  }

  std::cout << "\nl: ";
  for (int i=0; i<mN; i++) {
    std::cout << mLower[i] << " ";
  }

}

// Solves system of equations
void Spline::solveTridiaognal() {

  //Create delta and Gvec vectors of the Triangular system
  double *delta, *Gvec;
  delta = new double[mN+1];
  Gvec = new double[mN+1];
  for(int i=0; i<=mN; i++)
  {
  delta[i] = mDiag[i];
  Gvec[i] = mFvec[i];
  }

  // Elimination stage
  for(int i=1; i<=mN; i++)
  {
    delta[i] = delta[i] - mUpper[i-1]*(mLower[i-1]/delta[i-1]);
    Gvec[i] = Gvec[i] - Gvec[i-1]*(mLower[i-1]/delta[i-1]);

  }

  //Backsolve
  mCoeff[mN] = Gvec[mN]/delta[mN];
  mFullC[mN+1]=mCoeff[mN];
  for(int i=mN-1; i>=0; i--)
  {
    mCoeff[i] = ( Gvec[i] - mUpper[i]*mCoeff[i+1] )/delta[i];
    mFullC[i+1] = mCoeff[i];
  }
  mFullC[0] = mCoeff[1]-(double(1)/double(3))*mH*(*mFunction).derivative(mNodes[0]);
  mFullC[mN+2] =mCoeff[mN-1]+(double(1)/double(3))*mH*(*mFunction).derivative(mNodes[mN]);


  // Deallocates storage
  delete[] delta;
  delete[] Gvec;
}

// Alternative method for solving system of equations
void Spline::solveMethod2() {

  // Modified slightly to work for class

  char flag = 't';
  // Creates helper vectors to solve the system
  double* helperD = new double[mN+1];
  double* helperF = new double[mN+1];


  for(int i=0; i<=mN; i++)
  {
  helperD[i] = mDiag[i];
  helperF[i] = mFvec[i];
  }

  // elimination step
  for(int i=1; i<=mN; i++)
      {
          helperD[i] -= mUpper[i-1] * (mLower[i] / helperD[i-1]);
          helperF[i] -= helperF[i-1] * (mLower[i] / helperD[i-1]);
          if(helperD[i] == 0)
          {
              // Does not divide by 0
              flag = 'f';
              break;
          }

      }

  if(flag == 't')
  {
      // backsolve step
      mCoeff[mN] = helperF[mN]/helperD[mN];
      mFullC[mN+1]=mCoeff[mN];
      for(int i=mN-1; i>=0; i--)
          {

              mCoeff[i] = (helperF[i] - mUpper[i] * mCoeff[i+1]) / helperD[i];
              mFullC[i+1] = mCoeff[i];
          }
  }

  else
  {
      std::cout << "d = 0" << "\n";
  }

  mFullC[0] = mCoeff[1]-(double(1)/double(3))*mH*(*mFunction).derivative(mNodes[0]);
  mFullC[mN+2] =mCoeff[mN-1]+(double(1)/double(3))*mH*(*mFunction).derivative(mNodes[mN]);

  delete[] helperD;
  delete[] helperF;

}

// Shows coefficients of spline we have found
void Spline::showCoeff() {
  std::cout << "\nc: ";
  for (int i=0; i<mN+1; i++) {
    std::cout << mCoeff[i] << " ";
  }
}

// Evaluates the B function (which defines the cubic spline)
double Spline::evaluateB(const double xstar) {
  double B;
  if ((xstar <= -1)&&(xstar>=-2)) {
    B = pow((xstar+2), 3);
  } else if ((xstar >= 1)&&(xstar <=2)) {
    B = pow((2-xstar), 3);
  } else if ((xstar <= 0)&&(xstar>=-1)) {
    B = 1+(3*(xstar+1)) + (3*pow(xstar+1, 2)) - (3*pow(xstar+1,3));
  } else if ((xstar >= 0)&&(xstar <=1)) {
    B = 1+(3*(1-xstar)) + (3*pow(1-xstar, 2)) - (3*pow(1-xstar,3));
  } else {
    B = 0;
  }
  return B;
}

// Evaluate spline at point x
double Spline::evaluateSpline(const double x) {
  double sum = 0;
  double xstar;

  for(int i=0; i<=mN+2; i++) {
    xstar = (x-mFullNodes[i])/double(mH);
    sum = sum +(mFullC[i]*Spline::evaluateB(xstar));
  }

  return sum;
}

// Find absolute error of spline approximation to f at point x
double Spline::error(const double x) {
  double approx = Spline::evaluateSpline(x);
  double f = (*mFunction).evaluateF(x);
  return fabs(approx-f);
}
