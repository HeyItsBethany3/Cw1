#include "LSQ.hpp"
#include "AbstractFunction.hpp"
#include <cmath>
#include <iostream>
#include <cassert>
#include <string>
#include <fstream>


LSQ::LSQ(const int num, AbstractFunction& aFunction) {
  assert(num<=4); // n cannot be greater than 4
  bArray = new double[4];
  cArray = new double[num];
  cFull = new double[4];
  n = num;
  mFunction = &aFunction;
}

LSQ::~LSQ() {

  delete[] bArray;
  delete[] cArray;
  delete[] cFull;

  // TODO: delete mFunction?
}

double LSQ::phi_1(const double x)
{
  return(1.0);
}

double LSQ::phi_2(const double x)
{
  return(x);
}

double LSQ::phi_3(const double x)
{
  return(0.5*(3.0*pow(x,2)-1));
}

double LSQ::phi_4(const double x)
{
  return(0.5*(5.0*pow(x,3) - 3.0*x));
}


double LSQ::Gauss3(double (*P_function)(const double x))
{
  double sum = (double(8)/double(9))*(*mFunction).evaluateF(0)*(*P_function)(0);
  double xVal = sqrt(double(3)/double(5));
  sum += (double(5)/double(9))*((*mFunction).evaluateF(xVal)*(*P_function)(xVal)+(*mFunction).evaluateF(-xVal)*(*P_function)(-xVal));
  return sum;
}

double LSQ::Gauss5(double (*P_function)(const double x)) {
  double sum = (double(128)/double(225))*(*mFunction).evaluateF(0)*(*P_function)(0);
  double xVal = (double(1)/double(3))*sqrt(5.0-(2.0*sqrt(double(10)/double(7))));
  double w1 = double(322+13*double(sqrt(70)))/double(900);
  double xVal2 = (double(1)/double(3))*sqrt(5.0+(2.0*sqrt(double(10)/double(7))));
  double w2 = double(322-13*double(sqrt(70)))/double(900);
  sum += w1*((*mFunction).evaluateF(xVal)*(*P_function)(xVal)+(*mFunction).evaluateF(-xVal)*(*P_function)(-xVal));
  sum += w2*((*mFunction).evaluateF(xVal2)*(*P_function)(xVal2)+(*mFunction).evaluateF(-xVal2)*(*P_function)(-xVal2));
  return sum;
}

// Find B using gauss3 formula
void LSQ::findBGauss3() {
  bArray[0] = LSQ::Gauss3(LSQ::phi_1);
  bArray[1] = LSQ::Gauss3(LSQ::phi_2);
  bArray[2] = LSQ::Gauss3(LSQ::phi_3);
  bArray[3] = LSQ::Gauss3(LSQ::phi_4);
}
/* Note: you can call a class method within another class method but in order
to use a function pointer to another class method, that method must be static
*/

// Find B using gauss5 formula
void LSQ::findBGauss5() {
  bArray[0] = LSQ::Gauss5(LSQ::phi_1);
  bArray[1] = LSQ::Gauss5(LSQ::phi_2);
  bArray[2] = LSQ::Gauss5(LSQ::phi_3);
  bArray[3] = LSQ::Gauss5(LSQ::phi_4);
}



double LSQ::altGauss3(double (*P_function)(const double x)) {
  double x0, x1, x2;

  x0 = sqrt(3) /sqrt(5);
  x1 = 0.;
  x2 = - x0;

  double sum = (8./9.)*(*mFunction).evaluateF(x1)*(*P_function)(x1);
  sum += (5./9.)*(*mFunction).evaluateF(x0)*(*P_function)(x0);
  sum += (5./9.)*(*mFunction).evaluateF(x2)*(*P_function)(x2);
  return sum;

}

double LSQ::altGauss5(double (*P_function)(const double x)) {
  return 0; // Need to implement method
}

void LSQ::altFindBGauss3() {
  bArray[0] = LSQ::altGauss3(LSQ::phi_1);
  bArray[1] = LSQ::altGauss3(LSQ::phi_2);
  bArray[2] = LSQ::altGauss3(LSQ::phi_3);
  bArray[3] = LSQ::altGauss3(LSQ::phi_4);
}

void LSQ::altFindBGauss5() {
  bArray[0] = LSQ::altGauss5(LSQ::phi_1);
  bArray[1] = LSQ::altGauss5(LSQ::phi_2);
  bArray[2] = LSQ::altGauss5(LSQ::phi_3);
  bArray[3] = LSQ::altGauss5(LSQ::phi_4);
}




// Finds coefficients of the LSQ approx polynomial q_n
void LSQ::computeCoefficients() {
  double *Mu_inv_diag; //Vector holding the diagonal values of Mu inverse matrix

  Mu_inv_diag = new double[4];

  // Values were found analytically
  Mu_inv_diag[0] = 0.5;
  Mu_inv_diag[1] = 1.5;
  Mu_inv_diag[2] = 2.5;
  Mu_inv_diag[3] =double(7)/double(2);


  for(int i=0; i<4; i++) {
    cFull[i] = Mu_inv_diag[i] * bArray[i];
  }
  for(int i=0; i<n; i++)
  {
    cArray[i] = cFull[i];
  }

  delete[] Mu_inv_diag;
}

void LSQ::showB() {
  for(int i=0; i<n; i++)
  {
    std::cout << "b[" << i+1 << "] = " << bArray[i] << "\n";
  }
}

void LSQ::showC() {
  for(int i=0; i<n; i++)
  {
    std::cout << "c[" << i+1 << "] = " << cArray[i] << "\n";
  }
}

double LSQ::evaluateQ(const double x) {

  double* product = new double[4];
  product[0] = phi_1(x)*cFull[0];
  product[1] = phi_2(x)*cFull[1];
  product[2] = phi_3(x)*cFull[2];
  product[3] = phi_4(x)*cFull[3];
  //std::cout << "\nphi1: " << phi_1(x) << "\nphi2: "<< phi_2(x) << "\nphi3: " << phi_3(x) << "\nphi4: "<< phi_4(x) << "\n";


  double sum = 0;
  for(int i=0; i<n; i++) {
    sum += product[i];
  }
  //std::cout << "\nproduct1: " << product[0] << "\nproduct2: "<< product[1] << "\nproduct3: " << product[2]<< "\nproduct4: "<< product[3]<< "\n";
  delete[] product;
  return sum;
}

void LSQ::plotQ(const int nodes) {

  // Writes to file
  std::ofstream file;
  file.open("q.csv");
  assert(file.is_open());

  double x, q, exact, error;
  file << nodes << ","<<std::endl;
  double* fvec = new double[nodes+1];

  for (int i=0; i<=nodes; i++) {
    x =-1.0 +i*(2.0/double(nodes));
    q = LSQ::evaluateQ(x);
    fvec[i] = (*mFunction).evaluateF(x);
    file << q << ",";
  }
  file << std::endl;

  for (int i=0; i<=nodes; i++) {
    file << fvec[i] << ",";
  }
  delete[] fvec;

  file.close();
  std::system("rm ../../../MATLAB/q.csv");
  std::system("cp q.csv ../../../MATLAB");
}

double LSQ::errorNorm(const int nodes) {

  double x, q, exact, error;
  double norm = 0;

  for (int i=0; i<=nodes; i++) {
    x =-1.0 +i*(2.0/double(nodes));
    q = LSQ::evaluateQ(x);
    exact = (*mFunction).evaluateF(x);
    error = fabs(q-exact);
    norm += pow(error, 2);

  }
  return sqrt(norm);
}
