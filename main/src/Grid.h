#ifndef GRID
#define GRID

#include <vector>
#include <cmath>

using Pos = std::pair<int, int>;
using Vec2D = std::pair<double, double>;

template<typename Val>
class Grid {
  private:
    int size;
    std::vector< std::vector<Val> > grid;

  public:

    Grid(const int size, const Val init);
    Val              get(int i, int j);
    void             set(int i, int j, const Val &v);
    Val              operator()(int i, int j);
    void             operator()(int i, int j, const Val &v);
    bool             check_array_bound(int i, int j);
    std::vector<Pos> get_neighborhood(int i, int j);
    int              count_nn(const int i, const int j, const Val &val);
    double           curvature(int i, int j, const Val &occupied);
};


template<typename Val>
Grid<Val>::Grid(const int size, const Val init) :
  size(size),
  grid(std::vector< std::vector<Val> >(size, std::vector<Val>(size, init)))
{
}

template<typename Val>
Val Grid<Val>::get(int i, int j) {
  return this->grid[i][j];
}

template<typename Val>
void Grid<Val>::set(int i, int j, const Val &v) {
  this->grid[i][j] = v;
}


template<typename Val>
Val Grid<Val>::operator()(int i, int j) {
  return this->get(i, j);
}

template<typename Val>
void Grid<Val>::operator()(int i, int j, const Val &v) {
  this->set(i, j, v);
}


template<typename Val>
bool Grid<Val>::check_array_bound(int i, int j) {
  int size = this->size;
  return (0 <= i) && (i < size) && (0 <= j) && (j < size);
}

template<typename Val>
std::vector<Pos> Grid<Val>::get_neighborhood(int i, int j) {
  const bool even = i%2 == 0;

  /* hexagonal grid
   ___ ___ ___ ___ ___
  [___|___|___|___|___|
    [___|___|___|___|___|
  [___|___|___|___|___|
    [___|___|___|___|___|
  [___|___|___|___|___|

  */
  if (even) {
    return { Pos(i, j-1), Pos(i, j+1), Pos(i-1, j-1), Pos(i-1, j), Pos(i+1, j-1), Pos(i+1, j) };
  } else {
    return { Pos(i, j-1), Pos(i, j+1), Pos(i-1, j), Pos(i-1, j+1), Pos(i+1, j), Pos(i+1, j+1) };
  }
  // square grid
//  return { Pos(i-1, j), Pos(i+1, j), Pos(i, j-1), Pos(i, j+1) };
}

template<typename Val>
int Grid<Val>::count_nn(const int i, const int j, const Val &val) {
  const auto nn = this->get_neighborhood(i, j);
  int count = 0;

  for(auto& var : nn ) {
    if ( this->grid[var.first][var.second] == val) {
      count++;
    }
  }
  return count;
}

template<typename Val>
double Grid<Val>::curvature(int i, int j, const Val &occupied) {
  int count = this->count_nn(i, j, occupied);
  return (4-count)/3.0;
}


#endif
