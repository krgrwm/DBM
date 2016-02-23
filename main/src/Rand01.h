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
  private:
    std::default_random_engine generator;
    std::discrete_distribution<double> distribution;
  public:
    Pick(const std::vector<double>& vec) : distribution(vec.begin(), vec.end()) {}
    int operator()() {
      return this->distribution(generator);
    }
};

#endif
