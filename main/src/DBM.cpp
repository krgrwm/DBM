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
#include <chrono>

#include <stdio.h>

/* Private */

int DBM::center() {
  return this->b.center;
}

double DBM::calc_cluster_potential(const int i, const int j) {
  double sigma     = this->sigma;
  double curvature = this->b.cluster.curvature(i, j, true);
  double T_M       = this->b.val_cluster;

  return T_M * (1 - sigma*curvature);
}

void DBM::set_cluster_potential() {
  double new_value = 0.0;

  for (int i = 0; i < this->size; i++) {
    for (int j = 0; j < this->size; j++) {
      // !this->b.outer(i, j) は不要 (clusterとouterが排他的なら)
      // clusterが大きくなりすぎた場合を考えて一応条件に含めている
      if (this->b.cluster(i, j) && (!this->b.outer(i, j))) {
        new_value = calc_cluster_potential(i, j);
        this->grid(i, j, new_value);
      }
    }
  }
}

//void DBM::normalize_sigma() {
//  this->sigma = this->b.val_cluster / this->normalization_sigma * this->sigma;
//}

/* Public */
/* DEBUG */
//DBM::DBM(const int size, const double eta, const int N, const int threshold, const double sigma, SOR sor):
DBM::DBM(const int size, const double eta, const int N, const int threshold, const double sigma, SOR sor):
  size(size),
  eta(eta),
  N(N),
  sor(sor),
  r(),
  grid(size, 0.0), 
  __carvature(size, 0.0),
  b(size, 0.0, 100.0),
  peri(),
  threshold(threshold),
  counter(size, 0),
  sigma(sigma)
{
  std::random_device seedgen;
  this->mt = std::mt19937(seedgen());
  cout <<
    " size: "      << this->size               <<
    " N: "         << this->N                  <<
    " threshold: " << this->threshold          <<
    " sigma: "     << this->sigma              <<
    " eta: "       << this->eta                <<
    " omega: "     << this->sor.get_omega()    <<
    " epsilon: "   << this->sor.get_epsilon()  <<
    endl;
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
      add = !(this->b.is_boundary(i, j));
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
  const int c = this->b.center;
  add_particle(Pos(c, c));

  // set circular boundary (phi=1)
  for (int i = 0; i < this->size; i++) {
    for (int j = 0; j < this->size; j++) {
      if (this->b.is_outer_boundary(i, j)) {
        this->grid(i, j, this->b.val_outer);
        this->b.outer(i, j, true);
      }
    }
  }

  count = this->solve();
  return count;
}

int DBM::solve() {
//  this->sor._solve(this->size, this->grid, this->b); // for DEBUG
  int count;

  this->set_cluster_potential();
  count = this->sor.solve(this->size, this->grid, this->b);
  return count;
}


double DBM::grad_phi(const Pos& pos) {
  int i = pos.first;
  int j = pos.second;

//  /* DEBUG */
//  double grad = this->grid.grad_abs(i, j);
//  return grad;

  const auto nn = this->grid.get_neighborhood(i, j);

  // 最もpotentialの値が高い境界を探す
  // gradが最も大きくなるため
  double max=0.0;
  int _i, _j;
  for(auto& var : nn ) {
    _i = var.first;
    _j = var.second;
    if (this->b.cluster(_i, _j) && this->grid(_i, _j) > max) {
      max = this->grid(_i, _j);
    }
  }

//  return this->grid(i, j) - this->b.val_cluster;
  double grad_max = fabs(max - this->grid(i, j));

//  cout << "grad=" << grad << "  grad_max=" << grad_max << endl;
  return grad_max;
}

// calculate probability from potential
PList DBM::plist(Perimeter& peri) {
  cout << "plist" << endl;

//  double C = 0.0;  // Normalization constant

  PList plist(peri.size());

//  // calc normalization constant C
//  for (const auto& pos : peri) {
////    cout << "grad_phi=" << this->grid(pos.first, pos.second) << endl;
////    cout << "A*grad_phi=" << A * this->grid(pos.first, pos.second) << endl;
//    C += pow(A*this->grad_phi(pos), this->eta);
//  }
//  cout << "C=" << C << endl;
  // calc probability
  int i=0;
  for (const auto& pos : peri) {
//    cout << this->grad_phi(pos) << endl;
//    double p = pow(A*this->grad_phi(pos), this->eta) / C ;
    double p = pow(this->grad_phi(pos), this->eta);
    plist.append(i, pos, p);
    i++;
  }
  return plist;
}

PosVal DBM::select(PList& pl) {
  cout << "select" << endl;
  
  Pick pick(pl.get_plist());

//  /* DEBUG */
//  auto ps    = pl.get_plist();
//  auto picks = pick.distribution.probabilities();
//  std::sort(ps.begin(), ps.end());
//  std::sort(picks.begin(), picks.end());

//  int N=1000;
//  std::vector<int> n(ps.size());
//  for (int j=0; j<N; j++) {
//    n[pick(this->mt)] += 1;
//  }
//  cout << "DEBUG" << endl;
//  for (int j=0; j<ps.size(); j++) {
//    cout << picks[j] << ":" << double(n[j])/N << endl;
//  }


//  for (int j=0; j<ps.size(); j++) {
//    cout << ps[j] << ":" << picks[j] << endl;
//  }

  int i;
  Pos p;
  int pi, pj;
  int c;

  while (true) {
    i = pick(this->mt);
//    this->counter[pl.pos(i)] += 1;
    p = pl.pos(i);
    pi = p.first;
    pj = p.second;
    c = this->counter(pi, pj);
    this->counter(pi, pj, c+1);
    if (this->counter(pi, pj) >= this->threshold) {
      this->counter(pi, pj, 0);
      return pl.at(i);
    }
  }
}


//  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//  shuffle(pl.plist.begin(), pl.plist.end(), std::default_random_engine(seed));
//
//  double p;
//  double sum;
//  bool   count;
//
//  while (true) {
//    p     = this->r.rand();
//    sum   = 0;
//    count = false;
//
//    // sum >= p が満たされ、counterをincrementしたらfor文を抜け
//    // 再びpを設定しforを繰り返す
//    for (int i=0; i < pl.size() && !count; i++) {
//      sum += pl.p(i);
//      if (sum >= p) {
//        this->counter[pl.pos(i)] += 1;
//        if (this->counter[pl.pos(i)] >= this->threshold) {
//          // 0にする必要はない
//          // boundaryとなり候補に加わらなくなるため
//          this->counter[pl.pos(i)] = 0;
//          return pl.at(i);
//        }
//        count = true;
//      }
//    }
//  }
//}

void DBM::update_perimeters(const Pos& pos) {
  const auto new_peri = get_perimeter(pos);

//  for(auto& var : new_peri ) {
//    this->counter(var.first, var.second, 0);
//  }

  this->peri.insert(new_peri.begin(), new_peri.end());
}

void DBM::add_particle(const Pos& p) {
  cout << "add_particle" << endl;

  this->b.cluster(p.first, p.second, true);
  this->grid(p.first, p.second, this->b.val_cluster);
  update_perimeters(p);

  // delete site to stick from perimeters
  this->peri.erase(p);
}

void DBM::write_header(ofstream &ofs) {
  ofs <<
    "# size: "      << this->size                        <<
    " N: "          << this->N                           <<
    " eta: "        << double(this->eta)                 <<
    " threshold: "  << this->threshold                   <<
    " sigma: "      << double(this->sigma)               <<
    " omega: "      << double(this->sor.get_omega())     <<
    " epsilon: "    << double(this->sor.get_epsilon())   <<
    endl;
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
  write_header(gofs);
  write_header(bofs);
  write_header(hofs);

  // write data
  for (int i = 0; i < this->size; i++) {
    for (int j = 0; j < this->size; j++) {
      gofs << this->grid(i, j) << " ";
      bofs << this->b.cluster(i, j) << " ";

      if (this->b.cluster(i, j)) {
        newj = j + ((i%2)-0.5)/2.0;
        hofs << newj << " " << i << endl;
      }
    }
    gofs << endl;
    bofs << endl;
  }
}

Pos DBM::select_from_perimeters() {
  auto pl = this->plist(this->peri);
  const auto stick_pos = this->select(pl).first;
  return stick_pos;
}

void DBM::step() {
  this->add_particle(select_from_perimeters());

  cout << "solve" << endl;
  this->solve();
}
