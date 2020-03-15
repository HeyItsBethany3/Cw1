For Q2 part 1b, open the folder 'Q2b-individual'. Run the files individually.
The scripts contain option I ('Gaussian.cpp' and 'IoannaGaussian.cpp') and
option II ('LSQ_coeffs.cpp') in the form of functions, and they are tested.

For the rest of Q2 we have designed a class system. There is an 'AbstractFunction'
class, from which we have two derived classes: 'Function1' which implements the
function specified in the question and 'Function2' which implements a simple
testing function f(x)=1+x+x^2.

The important class is 'LSQ'. In order to construct a least squares approximation
p(x) you must pass in an 'AbstractFunction' object which defines f(x) and the order
of the approximation n. Then you construct the b=<f,p> coefficients using a
quadrature method, then compute the coefficients of the approximation. You can
evaluate p at a specific x value and find its error.
This is done in our 'Driver.cpp' file which is our main cpp file.

Compile and run all the files by using the makefile. Navigate to this directory
using the terminal and then type 'make'. Use 'make clean' to restore the directory
back to its original state.

We also create a plot of q(x) and f(x) in 'Driver.cpp', using a method 'plotQ'
of 'LSQ'. This creates 'q.csv'. We create the plot by running 'plotQ'.m in our
matlab working directory.
