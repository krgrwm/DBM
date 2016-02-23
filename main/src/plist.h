#include <vector>
#include <map>

using Pos       = std::pair<int, int>;
using PosVal    = std::pair<Pos, double>;

class PList {
//  private:
  public:
    std::vector<Pos>    poslist;
    std::vector<double> plist;
  public:
    PList(const int size);
    void append(const int i, const Pos &pos, const double p);
    int size() const;
    double p(const int i) const;
    Pos pos(const int i) const;
    PosVal at(const int i) const;
    const std::vector<double> get_plist();
};
