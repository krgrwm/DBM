#include "Grid.h"

class SOR {
  private:
    double omega;
    double epsilon;

  private:
    double sum(Grid<double> &grid, int i, int j);
  public:
    SOR();
    void solve(int size, Grid<double> &grid, Grid<bool> &boundary);
    void solve_max(int size, Grid<double> &grid, Grid<bool> &boundary);
};
