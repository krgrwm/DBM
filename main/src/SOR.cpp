#include "SOR.h"
#include <cmath>
#include <iostream>

SOR::SOR(const double omega, const double epsilon) : omega(omega), epsilon(epsilon){} 

int SOR::solve(int size, Grid<double> &grid, Boundary &boundary) {
  return this->_solve_max(size, grid, boundary);
}

//void SOR::_solve(int size, Grid<double> &grid, Grid<bool> &boundary) {
//  double gij       = 0.0;
//  double new_gij   = 0.0;
//  double sum       = 0.0;
//
//  double error     = 0.0;
//  const double epsilon = 1.0E-5;
//
//  do {
//    error = 0.0;
//    for (int i = 1; i < size-1; i++) {
//      for (int j = 1; j < size-1; j++) {
//        if ( !boundary(i, j)) {
//          gij = grid(i, j);
//          sum = this->sum(grid, i, j);
//          new_gij = gij + omega * ( sum /6.0 - gij );
//          grid(i, j, new_gij);
//
//          // calc error
//          error = error + fabs(new_gij - gij);
//        }
//      }
//    }
//    std::cout << "err " << error << ": " << (size-2)*(size-2)*epsilon << std::endl;
//    // errorの平均で判定
//  } while(error >= (size-2)*(size-2)*epsilon);
//  std::cout << "solve end" << std::endl;
//}

int SOR::_solve_max(int size, Grid<double> &grid, Boundary &boundary) {
  double gij       = 0.0;
  double new_gij   = 0.0;
  double sum       = 0.0;
  int count=0;

  double max       = 0.0;
  double err       = 0.0;

  do {
    max = 0.0;
    for (int i = 1; i < size-1; i++) {
      for (int j = 1; j < size-1; j++) {
        if ( !boundary.is_boundary(i, j) ) {
          gij = grid(i, j);
          sum = this->sum(grid, i, j);
          new_gij = gij + omega * ( sum /6.0 - gij );
          grid(i, j, new_gij);

          // calc error
//          err = fabs((new_gij - gij)/new_gij);
          err = fabs(new_gij - gij);
          max = fmax(err, max);
        }
      }
    }
    count++;
//    std::cout << "max " << max << ": " << this->epsilon << std::endl;
  } while(max >= this->epsilon);
  return count;
//  std::cout << "solve end" << std::endl;
}

double SOR::sum(Grid<double> &grid, int i, int j) {
  bool   even  = (i%2 == 0);

  if (even) {
    return grid(i, j-1) + grid(i, j+1) + grid(i-1, j-1) + grid(i-1, j) + grid(i+1, j-1) + grid(i+1, j);
  } else {
    return grid(i, j-1) + grid(i, j+1) + grid(i-1, j) + grid(i-1, j+1) + grid(i+1, j) + grid(i+1, j+1);
  }
}

double SOR::get_omega() const {
  return this->omega;
}

double SOR::get_epsilon() const {
  return this->epsilon;
}
