#ifndef ABSTRACTMATRIXHEADERDEF
#define ABSTRACTMATRIXHEADERDEF

class AbstractMatrix
{
  public:

    // Destructor
    virtual ~AbstractMatrix() = 0;

    //Matrix-vector multiplication
    double* operator*(const double* v) const;

    // Entries of matrix
    double mEntries[2][2];
    //TODO note: made public so could access matrix entries output from evaluateInvJacobi - ok? or make protected?


};



#endif
