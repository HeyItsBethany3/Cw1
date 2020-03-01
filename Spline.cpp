#include "Spline.hpp"
#include "AbstractFunction.hpp"

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

}



double Spline::GetLength() {
  return mFvec[0];
}
