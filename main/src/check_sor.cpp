#include <time.h>
#include <iostream>
#include <string>
#include "DBM.h"
#include <sys/types.h>
#include <unistd.h>
#include <map>


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
  const double eta       = atof(argv[3]);
  const int    threshold = atoi(argv[4]);
  const double sigma     = atof(argv[5]);
  const string file      = string(argv[6]);

  const double omega = 1.97;
  const double epsilon = 1E-3;
  auto dbm = DBM(size, eta, N, threshold, sigma, SOR(omega, epsilon));

  cout << "Initialize DBM" << endl;

  const int center = dbm.center();

  int r = 30;
  for (int i=center-r; i<center+r; i++) {
    for (int j=center-r; j<center+r; j++) {
      if (pow(i-center, 2) + pow(j-center, 2) < r) {
        dbm.add_particle(Pos(i, j));
      }
    }
  }

  int count = dbm.init();
  cout << count << endl;

  dbm.write(file);

  return 0;
}
