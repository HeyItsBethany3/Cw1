#include <iostream>
#include <cmath>

void euler(const double theta0, const double g, const double l, const int n) {
  const double alpha = sqrt(g/l);
  const double T = 7.4162987092/sqrt(alpha); // final time step
  const double h = T/double(n); // step-size
  double y1old = theta0;
  double y2old = 0;
  double y1, y2;
  for (int i=1; i<=n; i++) {
    y1 = y1old + (h*y2old);
    y2 = y2old + (-h*pow(alpha,2)*sin(y1old));
    y1old = y1;
    y2old = y2;
  }
  std::cout << "\ny1: " << y1 << " y2: " << y2 << "\n";


}

int main(int argc, char* argv[]) {
  const double theta0 = M_PI/double(2.0); // initial theta value
  euler(theta0, 1, 1, 100);

  return 0;
}
