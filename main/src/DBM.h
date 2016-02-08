#include "Rand01.h"
#include "Grid.h"
#include "plist.h"
#include "SOR.h"
#include "Boundary.h"

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
    bool           is_outer_interface(const int i, const int j);
    double         gibbs_thomson(const int i, const int j); // calculate change of potential (T)
    double         calc_cluster_potential(const int i, const int j); // calculate phi = phi0 - sigma/R
    void           set_cluster_potential();



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
    Boundary      b;
    Stick       stick; // stuck particle set
    Perimeter   peri;  // candidates to stick

    // noise-reduction
    int   threshold;
    map<Pos, int> counter;

    // surface tension
    double      sigma;
};

