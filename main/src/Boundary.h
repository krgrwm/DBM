#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Grid.h"

/* Circular Boundary */

class Boundary {
  private:

  public:
    Boundary(const int size, const double val_outer, const double val_cluster);
    bool is_boundary(const int i, const int j);
    bool is_outer_boundary(const int i, const int j);

  public:
    int size;
    const double val_outer;
    const double val_cluster;
    Grid<bool> outer;
    Grid<bool> cluster;
    int center; // position of seed
    int R_outer; // radius of outer circular boundary
};

#endif
