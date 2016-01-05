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


DBM::DBM(const int size, const double eta, const int N)
{
  this->eta   = eta;
  this->size  = size;
  this->grid  = Grid(size, vector<double>(size, 0.0));
  this->b     = Boundary(size, vector<bool>(size, false));
  this->stick = Stick();
  this->peri  = Perimeter();
  this->r     = Rand01();
  this->N     = N;
}

const bool DBM::check_bound(Pos p) {
  int i = p.first;
  int j = p.second;
  int size = this->size;
  return (0 <= i) && (i < size) && (0 <= j) && (j < size);
}

Perimeter DBM::get_perimeter(Pos p) {
  int i = p.first;
  int j = p.second;

  // neumann neighborhood
  const vector<Pos> candidates = { Pos(i-1, j), Pos(i+1, j), Pos(i, j-1), Pos(i, j+1) };
  Perimeter peris;

  // if pos is not boundary and inside grid -> add site as candidate
  auto f = [&](Perimeter& acc, const Pos& pos) 
  {
    const int i = pos.first;
    const int j = pos.second;
    const bool ok = check_bound(pos);
    bool add = false;

    if (ok) {
      add = !(this->b[i][j]);
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
  this->grid[c][c] = 0.0;
  this->b[c][c] = true;
  (this->stick).insert(Pos(c, c));

  // set circular boundary (phi=1)
  for (int i = 0; i < this->size; i++) {
    for (int j = 0; j < this->size; j++) {
      r2 = pow(j-c, 2) + pow(i-c, 2);
      if (r2 >= r*r) {
        this->grid[i][j] = 1.0;
        this->b[i][j] = true;
      }
    }
  }
  // add perimeters at (c, c)
  auto peris = get_perimeter(Pos(c, c));
  (this->peri).insert(peris.begin(), peris.end());

//  clock_t begin = clock();
  this->solve(1400);
//  clock_t end = clock();
//  cout << double(end-begin)/CLOCKS_PER_SEC << endl;
}

// SOR method
void DBM::solve(int N) {
  const auto omega = 1.5;
  double gij = 0.0;

  for (int n = 0; n < N; n++) {
    for (int i = 1; i < this->size-1; i++) {
      for (int j = 1; j < this->size-1; j++) {
        if ( !this->b[i][j] ) {
          gij = this->grid[i][j];
          this->grid[i][j] = gij + omega * ( 
              (this->grid[i+1][j] + this->grid[i-1][j] + this->grid[i][j+1] + this->grid[i][j-1])/4.0 - gij
              );
        }
      }
    }
  }
}

// calculate probability from potencial
vector<PosVal> DBM::plist(Perimeter& peri) {
  double C = 0.0;
  vector<PosVal> plist(peri.size());
  // calc normalization constant C
  for (const auto& pos : peri) {
    C += pow(this->grid[pos.first][pos.second], eta);
  }
  // calc probability
  int i=0;
  for (const auto& pos : peri) {
    double p = pow(this->grid[pos.first][pos.second], eta) / C ;
    plist[i] = PosVal(pos, p);
    i++;
  }
  return plist;
}

PosVal DBM::select(const vector<PosVal>& pl) {
  const double p = this->r.rand();
  double sum=0;

  for (int i=0; i < pl.size(); i++) {
    sum += pl[i].second;
    if (sum >= p) {
      return pl[i];
    }
  }
  return pl[pl.size()-1];
}

void DBM::update_perimeters(const Pos& pos) {
  const auto new_peri = get_perimeter(pos);
  this->peri.insert(new_peri.begin(), new_peri.end());
}

void DBM::add_particle(const PosVal& pv) {
  const Pos& p = pv.first;

  this->b[p.first][p.second] = true;
  this->grid[p.first][p.second] = 0;
  this->stick.insert(p);
  update_perimeters(p);
}

void DBM::write(const string& f) {
  const string gridfile = f + ".grid";
  const string boundaryfile = f + ".boundary";
  ofstream gofs(gridfile);
  ofstream bofs(boundaryfile);

  // write header
  gofs << "# size:" << this->size << " N:" << this->N << " eta:" << this->eta << endl;
  bofs << "# size:" << this->size << " N:" << this->N << " eta:" << this->eta << endl;

  // write grid data
  for(auto& row : this->grid ) {
    for(auto& var : row ) {
      gofs << var << " ";
    }
    gofs << endl;
  }

// write boudary data
  for(const auto& row : this->b) {
    for(const auto& var : row ) {
      bofs << var << " ";
    }
    bofs << endl;
  }
}

void DBM::step() {
  this->add_particle(this->select(this->plist(this->peri)));
  this->solve(50);
}

void DBM::grow() {
  for (int i=0; i<this->N; i++ ) {
    this->step();
  }
}
