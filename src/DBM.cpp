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


DBM::DBM(const int size, const double eta, const int N, const int threshold): grid(size, 0.0), b(size, false)
{
  this->eta       = eta;
  this->size      = size;
  this->stick     = Stick();
  this->peri      = Perimeter();
  this->r         = Rand01();
  this->N         = N;
  this->threshold = threshold;
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

void DBM::init()
{
  const auto c = int(this->size/2);
  const auto r = c-2;
  double r2=0;

  // set seed at center (phi=0)
  this->grid(c, c, 0.0);
  this->b(c, c, true);
  (this->stick).insert(Pos(c, c));

  // set circular boundary (phi=1)
  for (int i = 0; i < this->size; i++) {
    for (int j = 0; j < this->size; j++) {
      r2 = pow(j-c, 2) + pow(i-c, 2);
      if (r2 >= r*r) {
        this->grid(i, j, 1.0);
        this->b(i, j, true);
      }
    }
  }
  // add perimeters at (c, c)
  auto peris = get_perimeter(Pos(c, c));
  (this->peri).insert(peris.begin(), peris.end());

  this->solve();
}

// SOR method
void DBM::solve() {
  const auto omega = 1.5;
  double gij = 0.0;
  double new_gij = 0.0;
  double sum = 0.0;

  double error = 0.0;
  double error_sum = 0.0;
  double epsilon = 1.0E-5;

  do {
    error = 0.0;
    error_sum = 0.0;
    for (int i = 1; i < this->size-1; i++) {
      for (int j = 1; j < this->size-1; j++) {
        if ( !this->b(i, j)) {
          gij = this->grid(i, j);
          sum = 0.0;
          for(auto& var : this->grid.get_neighborhood(i, j) ) {
            sum += this->grid(var.first, var.second);
          }
          new_gij = gij + omega * ( sum /6.0 - gij );
          this->grid(i, j, new_gij);

          // calc error
          error = error + abs(new_gij - gij);
          error_sum = error_sum + abs(new_gij);
        }
      }
    }
//    cout << error/error_sum << endl;
  } while(error/error_sum >= epsilon);
}

// calculate probability from potential
PList DBM::plist(Perimeter& peri) {
  double C = 0.0; // Normalization constant
  PList plist(peri.size());

  // calc normalization constant C
  for (const auto& pos : peri) {
    C += pow(this->grid(pos.first, pos.second), eta);
  }
  // calc probability
  int i=0;
  for (const auto& pos : peri) {
    double p = pow(this->grid(pos.first, pos.second), eta) / C ;
    plist.append(i, PosVal(pos, p));
    i++;
  }
  return plist;
}

PosVal DBM::select(const PList& pl) {
  while (true) {
    double p = this->r.rand();
    double sum = 0;
    bool count=false;

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
  this->grid(p.first, p.second, 0.0);
  this->stick.insert(p);
  update_perimeters(p);

  // delete site to stick from perimeters
  this->peri.erase(p);
}

void DBM::write(const string& f) {
  const string gridfile = f + ".grid";
  const string boundaryfile = f + ".boundary";
  ofstream gofs(gridfile);
  ofstream bofs(boundaryfile);

  // write header
  gofs << "# size:" << this->size << " N:" << this->N << " eta:" << this->eta << endl;
  bofs << "# size:" << this->size << " N:" << this->N << " eta:" << this->eta << endl;

  // write data
  for (int i = 0; i < this->size; i++) {
    for (int j = 0; j < this->size; j++) {
      gofs << this->grid(i, j) << " ";
      bofs << this->b(i, j) << " ";
    }
    gofs << endl;
    bofs << endl;
  }
}

void DBM::write_hex(const string& f) {
  const string hexfile = f + ".hex";
  ofstream hofs(hexfile);
  double newj = 0;

  // write header
  hofs << "# size:" << this->size << " N:" << this->N << " eta:" << this->eta << endl;

  // write data
  for (int i = 0; i < this->size; i++) {
    for (int j = 0; j < this->size; j++) {
      if (this->b(i, j)) {
        newj = j + ((i%2)-0.5)/2.0;
          hofs << newj << " " << i << endl;
      }
    }
  }
}

void DBM::step() {
  this->add_particle(this->select(this->plist(this->peri)));
  this->solve();
}
