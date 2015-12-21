#include <ctime>
#include <iostream>
#include "DBM.h"

using namespace std;

int main(int argc, char const* argv[])
{
  // get parameter
  if (argc < 2) {
    cout << "N" << endl;
    return 0;
  }
  int N = atoi(argv[1]);

  auto dbm = DBM(300);

  cout << "Initialize DBM" << endl;
  dbm.init();

  cout << "0 - " << N << endl;
  for (int i=0; i<N; i++ ) {
    dbm.step();
  }

  dbm.write("/tmp/data");
  return 0;
}
