#include "Rand01.h"
#include <set>

using namespace std;

using Pos       = pair<int, int>;
using PosVal    = pair<Pos, double>;
using Grid      = vector< vector<double> >;
using Boundary  = vector< vector<bool> >;
using Stick     = set<Pos>;
using Perimeter = set<Pos>;

class DBM {
  private:
    Perimeter      get_perimeter(Pos p);
    const bool     check_bound(Pos p);
    vector<PosVal> plist(Perimeter& peri);
    PosVal         select(const vector<PosVal>& ps);
    void           update_perimeters(const Pos& pos);
    void           add_particle(const PosVal& pv);

  public:
    DBM(int size);
    void  init();
    void  solve(int n);
    void  write(const string& f);
    void  step();

  public:
    int         size;
    Grid        grid;
    Boundary    b;
    Stick       stick;
    Perimeter   peri;
    Rand01      r;

};
