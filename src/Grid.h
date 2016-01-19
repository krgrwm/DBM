#include <vector>

using Pos = pair<int, int>;

template<typename Val>
class Grid {
  private:
    int size;
    std::vector< std::vector<Val> > grid;
  public:

    Grid(const int size, const Val init);
    Val operator()(int i, int j);
    void operator()(int i, int j, const Val &v);
    bool check_array_bound(int i, int j);
    const vector<Pos> get_neighborhood(int i, int j);
};


template<typename Val>
Grid<Val>::Grid(const int size, const Val init) {
  this->size = size;
  this->grid = std::vector< std::vector<Val> >(size, std::vector<Val>(size, init));
}

template<typename Val>
Val Grid<Val>::operator()(int i, int j) {
  return this->grid[i][j];
}

template<typename Val>
void Grid<Val>::operator()(int i, int j, const Val &v) {
  this->grid[i][j] = v;
}

template<typename Val>
bool Grid<Val>::check_array_bound(int i, int j) {
  int size = this->size;
  return (0 <= i) && (i < size) && (0 <= j) && (j < size);
}

template<typename Val>
const vector<Pos> Grid<Val>::get_neighborhood(int i, int j) {
  bool even = i%2 == 0;

  // hexagonal grid
  if (even) {
    return { Pos(i, j-1), Pos(i, j+1), Pos(i-1, j-1), Pos(i-1, j), Pos(i+1, j-1), Pos(i+1, j) };
  } else {
    return { Pos(i, j-1), Pos(i, j+1), Pos(i-1, j), Pos(i-1, j+1), Pos(i+1, j), Pos(i+1, j+1) };
  }
//  return { Pos(i-1, j), Pos(i+1, j), Pos(i, j-1), Pos(i, j+1) };
}
