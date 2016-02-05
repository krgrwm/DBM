#include <iostream>
#include <string>
#include "DBM.h"
#include <sys/types.h>
#include <unistd.h>


using namespace std;

int main(int argc, char const* argv[])
{
  // read parameter
  if (argc < 6) {
    cout << "SystemSize N eta sigma output" << endl;
    return 0;
  }
  const int    size      = atoi(argv[1]);
  const int    N         = atoi(argv[2]);
  const double eta       = atoi(argv[3]);
  const int    threshold = atoi(argv[4]);
  const double sigma     = atoi(argv[5]);
  const string file      = string(argv[6]);


  double omega0 = 1.00;
  double omega = 0.0;
  double epsilon = 1E-5;
  int count = 0;

  cout << "omega iter" << endl;

  /* 1.00 to 1.70 */
  omega0 = 1.00;
  for (int i=0; i<8; i++) {
    omega = omega0 + i * 0.1;
    auto sor = SOR(omega, epsilon);

    cout << omega << " ";
    auto dbm = DBM(size, eta, N, threshold, sigma, sor);
    count = dbm.init();
    cout << count << endl;
}

/* 1.80 to 1.99 */
omega0 = 1.80;
for (int i=0; i<20; i++) {
  omega = omega0 + i * 0.01;
  auto sor = SOR(omega, epsilon);

  cout << omega << " ";
  auto dbm = DBM(size, eta, N, threshold, sigma, sor);
  count = dbm.init();
  cout << count << endl;
}


  return 0;
}
