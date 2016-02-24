#ifndef RAND01
#define RAND01

#include <random>
#include <vector>

using namespace std;

class Rand01 {
  private:
    mt19937                       mt;
    uniform_real_distribution<>   rand01;
  public:
    Rand01() {
      random_device rnd;
      this->mt     = mt19937(rnd());
      this->rand01 = uniform_real_distribution<>(0.0, 1.0);
    }
    double rand() {
      return this->rand01(mt);
    }
};

class RandInt {
  private:
    mt19937                         mt;
  public:
    RandInt(int N) {
      std::random_device rnddev;
      this->mt = mt19937(rnddev());
    }
    double rand(int start, int end) {
      std::uniform_int_distribution<int> rnd(start, end);
      return rnd(mt);
    }
};

class Pick {
//  private:
public:
    std::discrete_distribution<double> distribution;
  public:
    Pick(const std::vector<double>& vec) : distribution(vec.begin(), vec.end()) {}
    int operator()(std::mt19937 &mt) {
      return this->distribution(mt);
    }
};
class Pick_vicsec {
//  private:
public:
    std::vector<double> plist;
    std::uniform_int_distribution<int> disti;
    std::uniform_real_distribution<double> distr;
  public:
    Pick_vicsec(const std::vector<double>& plist) : plist(plist), disti(0, plist.size()-1) {
      double max = *std::max_element(plist.begin(), plist.end());
      distr = std::uniform_real_distribution<double>(0, max);
    }
    int operator()(std::mt19937 &mt) {
      int i;
      double rval;

      while (true) {
        i    = this->disti(mt);
        rval = this->distr(mt);

        if (rval < this->plist[i]) {
          return i;
        }
      }
    }
};

#endif
