#include "DBM.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include <numeric>
#include <algorithm>

/* Private */

bool DBM::is_outer_interface(const int i, const int j) {
  double r2 = pow(j-this->center, 2) + pow(i-this->center, 2);
  return (r2 >= this->R_circular_interface*this->R_circular_interface);
}

double DBM::gibbs_thomson(const int i, const int j) {
  double curvature = this->b.curvature(i, j, true);
  return this->sigma * curvature;
}

/* Public */
DBM::DBM(const int size, const double eta, const int N, const int threshold, const double sigma, SOR sor):
  size(size),
  eta(eta),
  N(N),
  sor(sor),
  r(),
  grid(size, 0.0), 
  b(size, false), 
  stick(),
  peri(),
  threshold(threshold),
  sigma(sigma),
  outer(0),
  cluster(100)
{
  this->center = int(size/2);
  this->R_circular_interface = this->center - 2;
}

Perimeter DBM::get_perimeter(Pos p) {
  int i = p.first;
  int j = p.second;

  const vector<Pos> candidates = this->grid.get_neighborhood(i, j);
  Perimeter peris;

  // if pos is not boundary and inside grid -> add site as candidate
  auto f = [&](Perimeter& acc, const Pos& pos) 
  {
    const int i = pos.first;
    const int j = pos.second;
    const bool ok = this->grid.check_array_bound(i, j);
    bool add = false;

    if (ok) {
      add = !(this->b(i, j));
    }

    if (add) {
      acc.insert(pos);
    }
    return acc;
  };

  peris = accumulate(candidates.begin(), candidates.end(), peris, f);
  return peris;
}

int DBM::init()
{
  int count=0;

  // set seed at center (phi=0)
  const int c = this->center;
  this->grid(c, c, 0.0);
  this->b(c, c, true);
  (this->stick).insert(Pos(c, c));

  // set circular boundary (phi=1)
  for (int i = 0; i < this->size; i++) {
    for (int j = 0; j < this->size; j++) {
      if (this->is_outer_interface(i, j)) {
        this->grid(i, j, this->outer);
        this->b(i, j, true);
      }
    }
  }
  // add perimeters at (c, c)
  auto peris = get_perimeter(Pos(c, c));
  (this->peri).insert(peris.begin(), peris.end());

  count = this->solve();
  return count;
}

int DBM::solve() {
//  this->sor._solve(this->size, this->grid, this->b); // for DEBUG
  return this->sor.solve(this->size, this->grid, this->b); // for DEBUG
}


double DBM::grad_phi(const Pos& pos) {
  int i = pos.first;
  int j = pos.second;
  double curvature = this->b.curvature(i, j, true);

//  return this->grid(i, j) + this->sigma*curvature;
  return this->grid(i, j) - this->cluster;
}

// calculate probability from potential
PList DBM::plist(Perimeter& peri) {
  double C = 0.0;  // Normalization constant

  // grad_phi << 1 -> C << 1 -> phi/C -> NaN
  // grad_pihが小さすぎる(1E-15)とかになるとCも小さくなり
  // phi/C = 0/0のようになるので定数Aをかけ大きな値にする
  // 特にeta=2以降だとこれが顕著になる
  double A = 1E15;
  PList plist(peri.size());

  // calc normalization constant C
  for (const auto& pos : peri) {
//    cout << "grad_phi=" << this->grid(pos.first, pos.second) << endl;
//    cout << "A*grad_phi=" << A * this->grid(pos.first, pos.second) << endl;
    C += pow(A*this->grad_phi(pos), eta);
  }
//  cout << "C=" << C << endl;
  // calc probability
  int i=0;
  for (const auto& pos : peri) {
    double p = pow(A*this->grad_phi(pos), eta) / C ;
    plist.append(i, PosVal(pos, p));
    i++;
  }
  return plist;
}

PosVal DBM::select(const PList& pl) {
//  cout << "select" << endl;
  while (true) {
    double p     = this->r.rand();
    double sum   = 0;
    bool   count = false;


    for (int i=0; i < pl.size() && !count; i++) {
      sum += pl.p(i);
      if (sum >= p) {
        this->counter[pl.pos(i)] += 1;
        if (this->counter[pl.pos(i)] >= this->threshold) {
          // 0にする必要はない
          // boundaryとなり候補に加わらなくなるため
          this->counter[pl.pos(i)] = 0;
          return pl.at(i);
        }
        count = true;
      }
    }
  }
}

void DBM::update_perimeters(const Pos& pos) {
  const auto new_peri = get_perimeter(pos);
  this->peri.insert(new_peri.begin(), new_peri.end());
}

void DBM::add_particle(const PosVal& pv) {
  const Pos& p = pv.first;

  this->b(p.first, p.second, true);
  this->grid(p.first, p.second, this->cluster);
  this->stick.insert(p);
  update_perimeters(p);

  // delete site to stick from perimeters
  this->peri.erase(p);
}

void DBM::write_header(ofstream &ofs) {
  ofs 
    << "# size:" << this->size 
    << " N:" << this->N 
    << " eta:" << this->eta 
    << " threshold:" << this->threshold
    << " sigma:" << this->sigma
    << " omega:" << this->sor.get_omega()
    << " epsilon:" << this->sor.get_epsilon() << endl;
}

void DBM::write(const string& f) {
  const string gridfile = f + ".grid";
  const string boundaryfile = f + ".boundary";
  const string hexfile = f + ".hex";
  double newj = 0;

  ofstream gofs(gridfile);
  ofstream bofs(boundaryfile);
  ofstream hofs(hexfile);

  // write header
  gofs << "# size:" << this->size << " N:" << this->N << " eta:" << this->eta << endl;
  bofs << "# size:" << this->size << " N:" << this->N << " eta:" << this->eta << endl;
  hofs << "# size:" << this->size << " N:" << this->N << " eta:" << this->eta << endl;

  // write data
  for (int i = 0; i < this->size; i++) {
    for (int j = 0; j < this->size; j++) {
      gofs << this->grid(i, j) << " ";
      bofs << this->b(i, j) << " ";

      if (this->b(i, j)) {
        newj = j + ((i%2)-0.5)/2.0;
          hofs << newj << " " << i << endl;
      }
    }
    gofs << endl;
    bofs << endl;
  }
}

void DBM::step() {
  cout << "add_particle" << endl;
  this->add_particle(this->select(this->plist(this->peri)));

  cout << "solve" << endl;
  this->solve();
}
