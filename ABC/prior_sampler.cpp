#include <string>
#include <iostream>
  using std::cout;
  using std::endl;
#include <iomanip>
  using std::setprecision;
#include <chrono>
#include <random>
#include <cstdlib>
  using std::atoi;
  using std::atof;
#include <ctime>
  using std::time;
#include <boost/random/mersenne_twister.hpp>
boost::mt19937 gent;

int main(int argc, char *argv[]){
  // Process arguments
  if (argc != 3)
  {
      std::cerr << "Usage: " << argv[0] << " Q_LOW Q_HIGH\n"
          "Sample from a log-uniform distribution with lower bound\n"
          "Q_LOW and upper bound Q_HIGH\n";
      return 1;
  }

  double q_low = atof(argv[1]);
  double q_high = atof(argv[2]);

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  gent.seed(seed);

  std::uniform_real_distribution<double> distribution(q_low, q_high);
  double raw_q_sampled = distribution(gent);

  double q_sampled = pow(10, raw_q_sampled);
  // Print sampled parameter
  std::cout.precision(17);
  std::cout << q_sampled << std::endl;

  return 0;
}
