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
  if (argc != 5)
  {
      std::cerr << "Usage: " << argv[0] << " Q_LOW Q_HIGH\n"
          "Sample from a log-uniform distribution with lower bound\n"
          "S_LOW and upper bound S_HIGH for the selection coefficient (s)\n"
          "Sample from a uniform distribution with lower bound\n"
          "H_LOW and upper bound H_HIGH for the dominance coefficient (h)\n";
      return 1;
  }

  double s_low = atof(argv[1]);
  double s_high = atof(argv[2]);

  double h_low = atof(argv[3]);
  double h_high = atof(argv[4]);

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  gent.seed(seed);

  std::uniform_real_distribution<double> s_distribution(s_low, s_high);
  double raw_s_sampled = s_distribution(gent);

  std::uniform_real_distribution<double> h_distribution(h_low, h_high);
  double h_sampled = h_distribution(gent);

  double s_sampled = pow(10, raw_s_sampled);
  // Print sampled parameter
  std::cout.precision(17);
  std::cout << s_sampled <<  " " << h_sampled << std::endl;

  return 0;
}
