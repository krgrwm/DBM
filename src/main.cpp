#include <time.h>
#include <iostream>
#include <string>
#include "DBM.h"
#include <sys/types.h>
#include <unistd.h>


using namespace std;

int main(int argc, char const* argv[])
{
  // read parameter
  if (argc < 5) {
    cout << "SystemSize N eta output" << endl;
    return 0;
  }
  const int    size = atoi(argv[1]);
  const int    N    = atoi(argv[2]);
  const double eta  = atoi(argv[3]);
  const string file = string(argv[4]);

  auto dbm = DBM(size, eta, N);

  cout << "Initialize DBM" << endl;
  dbm.init();

  cout << "0 - " << N << endl;

  for (int i=0; i<N; i++ ) {
    dbm.step();
    // print progress
    if (0 == i % (N/10)) {
      cout << getpid() << ": " << 100*i/N << endl;
    }
  }
  cout << getpid() << ": " << "DONE" << endl;
  dbm.write(file);

  return 0;
}
