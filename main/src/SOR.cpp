#include "SOR.h"
#include <cmath>
#include <iostream>

SOR::SOR(const double omega, const double epsilon) : omega(omega), epsilon(epsilon){} 

int SOR::solve(int size, Grid<double> &grid, Boundary &boundary, Grid<bool> &peri) {
  return this->_solve_max(size, grid, boundary, peri);
}

//int SOR::_solve_max(int size, Grid<double> &grid, Boundary &boundary) {
//  double gij       = 0.0;
//  double new_gij   = 0.0;
//  double sum       = 0.0;
//  int count=0;
//
//  double max       = 0.0;
//  double err       = 0.0;
//
//  do {
//    max = 0.0;
//    for (int i = 1; i < size-1; i++) {
//      for (int j = 1; j < size-1; j++) {
//        if ( !boundary.is_boundary(i, j) ) {
//          gij = grid(i, j);
//          sum = this->sum(grid, i, j);
//          new_gij = gij + omega * ( sum /4.0 - gij );
//          grid(i, j, new_gij);
//
//          err = fabs(new_gij - gij);
//          max = fmax(err, max);
//        }
//      }
//    }
//    count++;
//  } while(max >= this->epsilon);
//  return count;
//}

/* DEBUG */
int SOR::_solve_max(int size, Grid<double> &grid, Boundary &boundary, Grid<bool> &peri) {
  double gij       = 0.0;
  int count=0;

  double max       = 0.0;
  double err       = 0.0;

  do {
    max = 0.0;
    for (int i = 1; i < size-1; i++) {
      for (int j = 1; j < size-1; j++) {
        if (boundary.outer.grid[i][j]==false && boundary.cluster.grid[i][j]==false && peri.grid[i][j]==false) {
          gij = grid.grid[i][j];
          grid.grid[i][j] = 0.25 * (grid.grid[i+1][j]+grid.grid[i-1][j]+grid.grid[i][j+1]+grid.grid[i][j-1]);

          err = fabs(grid.grid[i][j] - gij);
          max = fmax(err, max);
        }
      }
    }
    count++;
  } while(max >= this->epsilon);
  return count;
}

double SOR::sum(Grid<double> &grid, int i, int j) {
  return grid(i-1, j) + grid(i+1, j) + grid(i, j-1) + grid(i, j+1);
}

double SOR::get_omega() const {
  return this->omega;
}

double SOR::get_epsilon() const {
  return this->epsilon;
}

