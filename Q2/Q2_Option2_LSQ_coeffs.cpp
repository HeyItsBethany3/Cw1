#include <iostream>
#include <cmath>

//TODO function prototypes

//Legendre Polynomials

double phi_1(const double x)
{
  return(1.0);
}

double phi_2(const double x)
{
  return(x);
}

double phi_3(const double x)
{
  return(0.5*(3.0*pow(x,2.0)-1));
}

double phi_4(const double x)
{
  return(0.5*(5.0*pow(x,3.0) - 3.0*x));
}




// LSQ coefficients system solver

double* solveLSQCoefficients(const int n, const double* b)
/*
Returns the coefficients of the LSQ approx polynomial q_n
Parameters: n = order of LSQ approx polynomial q_n
*/
{
  double *Mu_inv_diag; //Vector holding the diagonal values of Mu inverse matrix
  double *coeffs; //Vector of coefficients

  Mu_inv_diag = new double[4];
  coeffs = new double[4];

  // Values were found analytically
  Mu_inv_diag[0] = 0.5;
  Mu_inv_diag[1] = 1.5;
  Mu_inv_diag[2] = 2.5;
  Mu_inv_diag[3] = 9.0;

  for(int i=0; i<n; i++)
  {
    coeffs[i] = Mu_inv_diag[i] * b[i];
  }

  delete[] Mu_inv_diag;

  return coeffs;
}



// Testing Legendre Polynomials

void evaluateLegendre(double (*legendre)(double x))
{
  double point, output;

  for(int i=0; i<=8; i++)
  {
    point = -1+(i*0.25);
    output = (*legendre)(point);
    std::cout << "  at " << point << ": "<< output <<"\n";
  }
}

void testLegendre()
{
  std::cout << "P_1" <<"\n";
  evaluateLegendre(phi_1);

  std::cout << "\nP_2" <<"\n";
  evaluateLegendre(phi_2);

  std::cout << "\nP_3" <<"\n";
  evaluateLegendre(phi_3);

  std::cout << "\nP_4" <<"\n";
  evaluateLegendre(phi_4);
}


// int main(int argc, char* argv[])
// {
//   testLegendre();
// 
//   return 0;
// }
