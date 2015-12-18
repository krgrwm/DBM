#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <functional>
#include <numeric>
#include <algorithm>

using namespace std;
using Pos = pair<int, int>;
using PosVal = pair<Pos, double>;
using Grid = vector< vector<double> >;
using Boundary = vector< vector<bool> >;
using Stick = set<Pos>;
using Perimeter = set<Pos>;

class DBM {
  private:
    int size;
    Grid grid;
    Boundary b;
    Stick stick;
    Perimeter peri;

  private:
    Perimeter get_perimeter(Pos p);
    const bool check_bound(Pos p);
  public:
    DBM(int);
    void init();
    void solve(int n);
    vector<PosVal> plist(Perimeter& peri);
};

// TODO: perimeterをDBMに持たせ、stickに追加した時に新たなperimeterを追加するようにする

DBM::DBM(int size)
{
  this->grid  = Grid(size, vector<double>(size, 0.0));
  this->b     = Boundary(size, vector<bool>(size, false));
  this->stick = Stick();
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
  vector<Pos> candidates = { Pos(i-1, j), Pos(i+1, j), Pos(i, j-1), Pos(i, j+1) };
  Perimeter peris;

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
      double val = this->grid[i][j];
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

  this->grid[c][c] = 0.0;
  this->b[c][c] = true;
  (this->stick).insert(Pos(c, c));

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
}

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


vector<PosVal> DBM::plist(Perimeter& peri) {
  double eta=1.0;
  double C = 0.0;
  vector<PosVal> plist(peri.size());
  // ignore potential < 1.0e-10
  for (const auto& pos : peri) {
    C += this->grid[pos.first][pos.second]; // ここでは無視しない(値が小さいから無視できる)
  }
  for (const auto& pos : peri) {
    if (this->grid[pos.first][pos.second] > 1.0e-10) {
      // TODO
    }
  }
}

int main(int argc, char const* argv[])
{
  return 0;
}
