#include "Grid.h"

class SOR {
  private:
    double omega;
    double epsilon;

  private:
    double sum(Grid<double> &grid, int i, int j);
  public:
//    SOR();
    SOR(const double omega, const double epsilon);
    int solve(int size, Grid<double> &grid, Grid<bool> &boundary);
    int _solve_max(int size, Grid<double> &grid, Grid<bool> &boundary);
    double get_omega() const;
    double get_epsilon() const;
};
