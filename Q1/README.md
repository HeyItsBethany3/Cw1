For Q1 part 1b, open the folder 'Q1b-individual'. Run the files individually.
The scripts contain option I ('BrittaSolveTridiaognal.cpp' and 'IoannaSolveTridiagonal')
and option II ('completeSpline.cpp') in the form of functions, and they are tested.

For the rest of Q1 we have designed a class system. There is an 'AbstractFunction'
class, from which we have a derived class 'Function1', which implements the
function specified in the question. 'Function1' has no attributes but two methods
which compute f(x) and f'(x).
The important class is 'Spline'. In order to construct a spline you must pass in
an 'AbstractFunction' object which defines f(x). To find the spline, you must
create the interpolation nodes, find the system to solve and then solve the
tridiagonal system. This is done in our 'execute.cpp' file which is our main cpp
file.

Compile and run all the files by using the makefile. Navigate to this directory
using the terminal and then type 'make'. Use 'make clean' to restore the directory
back to its original state.

In our 'execute.cpp' file we create a table by outputting data to 'Table.csv'.
We also output data to 'PlotSpline.csv' which we use to create a plot, by running
'PlotSpline.m' in our MATLAB working directory.
