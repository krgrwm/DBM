#include "Boundary.h"
#include <cmath>

Boundary::Boundary(const int size, const double val_outer, const double val_cluster) : 
  size(size),
  val_outer(val_outer),
  val_cluster(val_cluster),
  outer(size, false), 
  cluster(size, false)
{
  this->center = int(size/2);
  this->R_outer = this->center - 2;

  // set seed
  this->cluster(center, center, true);
}

bool Boundary::is_boundary(const int i, const int j) {
  return this->outer(i, j) || this->cluster(i, j);
}

bool Boundary::is_outer_boundary(const int i, const int j) {
  double r2 = pow(j-this->center, 2) + pow(i-this->center, 2);
  return (r2 >= this->R_outer*this->R_outer);
}
