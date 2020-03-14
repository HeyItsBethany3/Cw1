#include "AbstractMatrix.hpp"

class TestJacobiInvMatrix: public AbstractMatrix
/*
Inverse Jacobi matrix from Lecture 3, slide 8.
*/
{

  public:
    // Specialised Constructor
    TestJacobiInvMatrix(const double* x);

    // Destructor
    ~TestJacobiInvMatrix();

};
