#include "Spline.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include<iostream>
#include<fstream>
#include<cmath>


int main(int argc, char* argv[]) {

  // Question 1c and 1d
  int n = 8;
  double len = 2.0; // So h=1/4

  // Britta's method (for solving the tridiagonal system)
  Function1* f1 = new Function1(); // Initialises function object
  Spline* s1 = new Spline(len,n,*f1); // Initialises spline object
  (*s1).Nodes(); // Creates nodes
  (*s1).FindSystem(); // Finds system of equations
  //(*s1).showSystem();
  (*s1).solveTridiaognal(); // Solves system of equations
  std::cout << "Method 1: ";
  (*s1).showCoeff(); // Displays spline coefficients

  // Ioanna's method (for solving the tridiagonal system)
  Function1* f2 = new Function1();
  Spline* s2 = new Spline(len,n,*f2);
  (*s2).Nodes();
  (*s2).FindSystem();
  (*s2).solveMethod2(); // Alternative method for solving system
  std::cout << "\nMethod 2: ";
  (*s2).showCoeff();

  // Checks spline is exact at interpolating points
  std::cout << "\n\nSpline at x=0: " << (*s1).evaluateSpline(0);
  std::cout << "\nf(0): " << (*f1).evaluateF(0) << std::endl;

  std::cout << "\nSpline at x=0.25: " << (*s1).evaluateSpline(0.25);
  std::cout << "\nf(0.25): " << (*f1).evaluateF(0.25) << std::endl;

  // Checks spline works for other points in interval [0,len]
  std::cout << "\nSpline at x=0.3: " << (*s1).evaluateSpline(0.3);
  std::cout << "\nf(0.3): " << (*f1).evaluateF(0.3) << std::endl;

  std::cout << "\nSpline at x=0.9: " << (*s1).evaluateSpline(0.9);
  std::cout << "\nf(0.9): " << (*f1).evaluateF(0.9) << std::endl;

  std::cout << "\nSpline at x=0.1: " << (*s1).evaluateSpline(0.1);
  std::cout << "\nf(0.1): " << (*f1).evaluateF(0.1) << std::endl;


  // Question 1e

  // Writes to file
  std::ofstream file; // Used to create table
  file.open("Table.csv");
  assert(file.is_open());

  std::ofstream file2; // Used to create plot
  file2.open("PlotSpline.csv");
  assert(file2.is_open());

  // Output to file (table format)
  file << "h," << "q," << "error" << std::endl;
  //file << std::scientific;

  double xstar = double(1)/double(3);
  int sim = 10; // Change this to get more results
  double* hValue = new double[sim];
  double* error = new double[sim];

  for(int i=1; i<=sim; i++) {
    hValue[i-1] = double(1)/double(pow(2,i));
    int n = len/hValue[i-1]; // only works with specific len (otherwise n is not an integer)
    Spline* s = new Spline(len,n,*f1);

    (*s).Nodes();
    (*s).FindSystem();
    (*s).solveTridiaognal();
    error[i-1] = (*s).error(xstar); // Finds error
    double approx = (*s).evaluateSpline(xstar); // Evaluates q (spline at x)

    file << hValue[i-1] << "," << approx << "," << error[i-1] << std::endl;
    delete s;
  }

  for (int i=1; i<=sim; i++) {
    file2 << hValue[i-1] << ","; // h (step-size)
  }
  file2 << std::endl;
  for (int i=1; i<=sim; i++) {
    file2 << error[i-1] << ","; // Error
  }
  // Close files
  file.close();
  file2.close();

  // Copies plot file to matlab directory so it is easy for matlab file to run
  std::system("rm ../../../MATLAB/PlotSpline.csv");
  std::system("cp PlotSpline.csv ../../../MATLAB");

  // Deallocates storage
  delete f1;
  delete s1;
  delete f2;
  delete s2;
  delete[] hValue;
  delete[] error;

  return 0;
}
