#include <random>

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
