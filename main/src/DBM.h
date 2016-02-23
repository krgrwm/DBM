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
using Perimeter = set<Pos>;

class DBM {
  private:
    Perimeter      get_perimeter(Pos p);
    PList          plist(Perimeter& peri);

    // select site to stick according to rule
    PosVal         select(PList& ps);

    // add new perimeters(candidates)
    void           update_perimeters(const Pos& pos);
    double         grad_phi(const Pos& pos);
    void           write_header(ofstream &ofs);
    bool           is_outer_interface(const int i, const int j);
    double         gibbs_thomson(const int i, const int j); // calculate change of potential (T)
    double         calc_cluster_potential(const int i, const int j); // calculate phi = phi0 - sigma/R



  public:
//    DBM(const int size, const double eta, const int N, const int threshold, const double sigma, SOR sor);
/* DEBUG */
    DBM(const int size, const double eta, const int N, const int threshold, const double sigma, SOR_Square sor);
    int init();
    int solve();
    void  write(const string& f);
    void  write_hex(const string& f);
    Pos   select_from_perimeters();
    void  step();

    int  center();
    void           set_cluster_potential();
    void           add_particle(const Pos& p);

  private:
    const int    size;
    const double eta;
    const int    N;     // steps
    /* DEBUG */
//    SOR          sor;
    SOR_Square          sor;

    Rand01      r;
//    Grid<double>  grid;
//    Grid<double>  __carvature;
/* DEBUG */
    Grid_Square<double>  grid;
    Grid_Square<double>  __carvature;
    Boundary      b;
    Perimeter   peri;  // candidates to stick

    // noise-reduction
    int   threshold;
    map<Pos, int> counter;

    // surface tension
    const double      sigma;

    // random engine
    std::mt19937 mt;
};

