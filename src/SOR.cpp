#include "SOR.h"
#include <cmath>
#include <iostream>

//SOR::SOR() : omega(1.25), epsilon(1.0E-5)
SOR::SOR() : omega(1.25), epsilon(1.0E-4)
{}

void SOR::solve(int size, Grid<double> &grid, Grid<bool> &boundary) {
  double gij       = 0.0;
  double new_gij   = 0.0;
  double sum       = 0.0;

  double error     = 0.0;
  double error_sum = 0.0;

  do {
    error = 0.0;
    error_sum = 0.0;
    for (int i = 1; i < size-1; i++) {
      for (int j = 1; j < size-1; j++) {
        if ( !boundary(i, j)) {
          gij = grid(i, j);
          sum = this->sum(grid, i, j);
          new_gij = gij + omega * ( sum /6.0 - gij );
          grid(i, j, new_gij);

          // calc error
          error = error + fabs(new_gij - gij);
          error_sum = error_sum + fabs(new_gij);
        }
      }
    }
//    cout << error/error_sum << endl;
  } while(error/error_sum >= this->epsilon);
}

void SOR::solve_max(int size, Grid<double> &grid, Grid<bool> &boundary) {
  double gij       = 0.0;
  double new_gij   = 0.0;
  double sum       = 0.0;

  double max       = 0.0;
  double err       = 0.0;

  do {
    max = 0.0;
    for (int i = 1; i < size-1; i++) {
      for (int j = 1; j < size-1; j++) {
        if ( !boundary(i, j)) {
          gij = grid(i, j);
          sum = this->sum(grid, i, j);
          new_gij = gij + omega * ( sum /6.0 - gij );
          grid(i, j, new_gij);

          // calc error
          err = fabs((new_gij - gij)/new_gij);
          max = fmax(err, max);
        }
      }
    }
//    std::cout << max << std::endl;
  } while(max >= this->epsilon);
}


double SOR::sum(Grid<double> &grid, int i, int j) {
  bool   even  = (i%2 == 0);

  if (even) {
    return grid(i, j-1) + grid(i, j+1) + grid(i-1, j-1) + grid(i-1, j) + grid(i+1, j-1) + grid(i+1, j);
  } else {
    return grid(i, j-1) + grid(i, j+1) + grid(i-1, j) + grid(i-1, j+1) + grid(i+1, j) + grid(i+1, j+1);
  }
}
