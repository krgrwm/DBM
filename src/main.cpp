#include <time.h>
#include <iostream>
#include <string>
#include "DBM.h"

using namespace std;

int main(int argc, char const* argv[])
{
  // get parameter
  if (argc < 5) {
    cout << "SystemSize N eta output" << endl;
    return 0;
  }
  const int size    = atoi(argv[1]);
  const int N       = atoi(argv[2]);
  const double eta  = atoi(argv[3]);
  const string file = string(argv[4]);

  auto dbm = DBM(size , eta, N);

  cout << "Initialize DBM" << endl;
  dbm.init();

  cout << "0 - " << N << endl;

  dbm.grow();

//  for (int i=0; i<N; i++ ) {
//    dbm.step();
//  }

  dbm.write(file);
  return 0;
}
