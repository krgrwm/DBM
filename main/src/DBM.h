#include "Rand01.h"
#include "Grid.h"
#include "plist.h"
#include "SOR.h"

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
    double         grad_phi(const Pos& pos);
    void           write_header(ofstream &ofs);



  public:
    DBM(const int size, const double eta, const int N, const int threshold, const double sigma, SOR sor);
    int init();
    int solve();
    void  write(const string& f);
    void  write_hex(const string& f);
    void  step();

  private:
    int         size;
    double      eta;
    int         N;     // steps
    SOR         sor;

    Rand01      r;
    Grid<double>  grid;
    Grid<bool>        b;
    Stick       stick; // stuck particle set
    Perimeter   peri;  // candidates to stick

    // noise-reduction
    int   threshold;
    map<Pos, int> counter;

    // surface tension
    double      sigma;
};

