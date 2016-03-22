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

  const int c = dbm.center();
  dbm.init();

  cout << c << endl;

  int size_cluster = 50;
  for (int i=0; i<size_cluster; i++) {
    for (int j=i-1; j<size_cluster-i; j++) {
      dbm.add_particle(Pos(c+i, c+j));
    }
  }
  dbm.set_cluster_potential();

//  int size_cluster = 20;
//  for (int i=c; i<c+size_cluster; i++) {
//    for (int j=c; j<c+size_cluster; j++) {
//      dbm.add_particle(Pos(i, j));
//    }
//  }

//  dbm.add_particle(Pos(c, c));
//  dbm.add_particle(Pos(c, c+1));
//  dbm.add_particle(Pos(c, c+2));
//  dbm.add_particle(Pos(c, c+3));
//  dbm.add_particle(Pos(c, c+4));
//  dbm.add_particle(Pos(c+1, c));
//  dbm.add_particle(Pos(c+1, c+1));
//  dbm.add_particle(Pos(c+1, c+2));
//  dbm.add_particle(Pos(c+1, c+3));
//  dbm.add_particle(Pos(c+2, c+1));
//  dbm.add_particle(Pos(c+2, c+2));
//  dbm.add_particle(Pos(c+2, c+3));
//  dbm.add_particle(Pos(c+3, c+1));
//  dbm.add_particle(Pos(c+3, c+2));

  for (int i=0; i<N; i++ ) {
    dbm.step();
    cout << getpid() << ": " << 100.0*double(i)/N << endl;
  }

  dbm.write(file);

  return 0;
}
