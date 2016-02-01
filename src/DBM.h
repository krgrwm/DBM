#include "Rand01.h"
#include "Grid.h"
#include "plist.h"

#include <set>
#include <map>

using namespace std;

using Pos       = pair<int, int>;
using PosVal    = pair<Pos, double>;
using Stick     = set<Pos>;
using Perimeter = set<Pos>;

class DBM {
  private:
    Perimeter      get_perimeter(Pos p);
    PList          plist(Perimeter& peri);

    // select site to stick according to rule
    PosVal         select(const PList& ps);

    // add new perimeters(candidates)
    void           update_perimeters(const Pos& pos);
    void           add_particle(const PosVal& pv);

  public:
    DBM(const int size, const double eta, const int N, const int threshold);
    void  init();
    void  solve(int n);
    void  write(const string& f);
    void  write_hex(const string& f);
    void  step();

  private:
    int         size;
    Grid<double>  grid;
    Grid<bool>        b;
    Stick       stick; // stuck particle set
    Perimeter   peri;  // candidates to stick
    Rand01      r;
    double      eta;
    int         N;     // steps

    // noise-reduction
    int   threshold;
    map<Pos, int> counter;
};

