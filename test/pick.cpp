#include <iostream>
#include "../main/src/Rand01.h"

using namespace std;

int main(int argc, char const* argv[])
{
  vector<double> v = {0.1, 0.3, 0.6};
  vector<int>    n = {0, 0, 0};
  Pick pick(v);

  int N=100000;

  for (int i=0; i<N; i++) {
    n[pick()] += 1;
  }

  cout << double(n[0])/N << endl;
  cout << double(n[1])/N << endl;
  cout << double(n[2])/N << endl;

  return 0;
}
