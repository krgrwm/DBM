#include <vector>
#include <map>

using Pos       = std::pair<int, int>;
using PosVal    = std::pair<Pos, double>;

class PList {
  private:
    std::vector<PosVal> plist;
  public:
    PList(const int size);
    void append(const int i, const PosVal pv);
    int size() const;
    double p(const int i) const;
    Pos pos(const int i) const;
    PosVal at(const int i) const;
};
