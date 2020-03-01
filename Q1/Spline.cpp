#include "Spline.hpp"
#include "AbstractFunction.hpp"
#include<iostream>

// Constructor
Spline::Spline(const double len, const int n, AbstractFunction& aFunction) {

  //TODO:: Checks on input
  mLen = len;
  mN = n;
  mNodes = new double[n+1];
  mH = len/double(n);

  mFvec = new double[n+1];
  mDiag = new double[n+1];
  mUpper = new double[n];
  mLower = new double[n];
  mCoeff = new double[n+1];

  mFunction = &aFunction;

}

// Destructor
Spline::~Spline() {
  delete mNodes;
  delete mFvec;
  delete mDiag;
  delete mUpper;
  delete mLower;
}

// Constructs interpolating nodes
void Spline::Nodes() {
  for (int i=0; i<mN+1; i++) {
    mNodes[i] = i*mH;
  }
}

// Finds system of equations for spline
void Spline::FindSystem() {
  for (int i=1; i<mN; i++) {
    mFvec[i]=(*mFunction).evaluateF(mNodes[i]);
  }
  mFvec[0]=(*mFunction).evaluateF(mNodes[0])+((double(1)/double(3))*mH*(*mFunction).derivative(mNodes[0]));
  mFvec[mN]=(*mFunction).evaluateF(mNodes[mN])-((double(1)/double(3))*mH*(*mFunction).derivative(mNodes[mN]));

  // Find diagonal elements
  for (int i=0; i<=mN; i++) {
    mDiag[i] = 4;
  }

  // Construct upper diagonal elements
  mUpper[0]=2;
  for (int i=1; i<mN; i++) {
    mUpper[i]=1;
  }

  // Construct lower diagonal elements
  mLower[mN-1]=2;
  for (int i=0; i<mN-1; i++) {
    mLower[i]=1;
  }

}

void Spline::solveTridiaognal() {
  // Elimination stage
  for(int i=1; i<=mN; i++)
  {
    mDiag[i] = mDiag[i] - mUpper[i-1]*(mLower[i-1]/mDiag[i-1]);
    mFvec[i] = mFvec[i] - mFvec[i-1]*(mLower[i-1]/mDiag[i-1]);
  }

  //Backsolve
  mCoeff[mN] = mFvec[mN]/mDiag[mN];
  for(int i=mN-1; i>=0; i--)
  {
    mCoeff[i] = ( mFvec[i] - mUpper[i]*mCoeff[i+1] )/mDiag[i];
  }
}


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

void Spline::showCoeff() {
  std::cout << "\nc: ";
  for (int i=0; i<mN+1; i++) {
    std::cout << mCoeff[i] << " ";
  }
}