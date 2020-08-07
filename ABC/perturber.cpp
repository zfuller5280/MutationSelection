#include <math.h>
#include <string>
#include <iostream>
#include <random>
#include <chrono>
//#include <boost/random.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/random/mersenne_twister.hpp>
//using namespace boost::math;
//double alpha, beta, randFromUnif;

using namespace boost::math;
boost::mt19937 gent;

double randnum (double a, double b)
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  gent.seed(seed);
  std::uniform_real_distribution<double> distribution (a,b);
  return distribution(gent);
}

int main(int argc, char *argv[])
{
  // Process arguments
  if (argc != 2)
  {
      std::cerr << "Usage: " << argv[0] << " SCALE\n"
          "Perturb given parameter by drawing "
          "from beta distribution\n"
          "with expected value equal to the given parameter "
          "and with scaling parameter SCALE\n";
      return 1;
  }

  double scale = std::stod(argv[1]);

  // Seed random number generator
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  //typedef boost::random::mt19937 RandomNumberGenerator;


  //RandomNumberGenerator Rng(seed);

  // Read generation t and parameter from stdin
  int t;
  std::cin >> t;

  double q;
  std::cin >> q;

  if (std::cin.fail())
  {
      std::cerr << "Error: could not read t or q from stdin\n";
      return 1;
  }

  // Perturb parameter
  std::cout.precision(17);
  double alpha = q * scale;
  double beta = (1-q) * scale;
  double randFromUnif = randnum(0, 1);
  beta_distribution<> d(alpha, beta);
  double q_perturbed = quantile(d, randFromUnif);

  // Print perturbed parameter
  std::cout << q_perturbed << std::endl;

  return 0;
}
