#include "Spline.hpp"
#include "AbstractFunction.hpp"
#include<iostream>
#include<cassert>

// Constructor
Spline::Spline(const double len, const int n, AbstractFunction& aFunction) {

  // Checks interval and n is reasonable
  assert(len>0);
  assert(n>0);

  mLen = len; // spline interval is [0,len]
  mN = n; // n+1 interpolating points
  mNodes = new double[n+1]; // stores interpolation points
  mH = len/double(n); // distance between nodes

  mFvec = new double[n+1]; // f values at each node
  mDiag = new double[n+1]; // diagonal vector
  mUpper = new double[n]; // upper diaagonal vector
  mLower = new double[n]; // lower diagonal vector
  mCoeff = new double[n+1]; // coefficients of spline

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
}

// Constructs interpolating nodes
void Spline::Nodes() {
  for (int i=0; i<mN+1; i++) {
    mNodes[i] = i*mH;
  }
}

// Finds system of equations for spline
void Spline::FindSystem() {
  // We are creating a system of the form Ac=b where A is a tridiagonal matrix,
  // b is mFvec and c is mCoeff

  // Find the RHS vector for the system of equations, which we call mFvec
  for (int i=1; i<mN; i++) {
    mFvec[i]=(*mFunction).evaluateF(mNodes[i]);
  }
  mFvec[0]=(*mFunction).evaluateF(mNodes[0])+((double(1)/double(3))*mH*(*mFunction).derivative(mNodes[0]));
  mFvec[mN]=(*mFunction).evaluateF(mNodes[mN])-((double(1)/double(3))*mH*(*mFunction).derivative(mNodes[mN]));

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
  for(int i=mN-1; i>=0; i--)
  {
    mCoeff[i] = ( Gvec[i] - mUpper[i]*mCoeff[i+1] )/delta[i];
  }

  // Deallocates storage
  delete[] delta;
  delete[] Gvec;
}

// Displays system to solve
void Spline::showSystem() {
  std::cout << "\nF(x): ";
  for (int i=0; i<=mN; i++) {
    std::cout << mFvec[i] << " ";
  }
  // TODO: Are my f values actually right?

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

// Shows coefficients of spline we have found
void Spline::showCoeff() {
  std::cout << "\nc: ";
  for (int i=0; i<mN+1; i++) {
    std::cout << mCoeff[i] << " ";
  }
}
