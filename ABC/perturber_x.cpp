#include <string>
#include <iostream>
#include <random>
#include <chrono>
#include <math.h>

int main(int argc, char *argv[])
{
    // Process arguments
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " STDEV1 STDEV2\n"
            "Perturb two parameters by drawing "
            "from log-normal distribution\n"
            "centered on the given parameter (s)"
            "and with standard deviation STDEV1\n"
            "and from noral distribution\n"
            "centered on the given parameter (h)"
            "and with standard deviation STDEV2\n";
        return 1;
    }

    double stdev1 = std::stod(argv[1]);
    double stdev2 = std::stod(argv[2]);

    // Seed random number generator
    unsigned seed =
        std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    // Read generation t and parameter from stdin
    int t;
    std::cin >> t;

    double s, h;
    std::cin >> s >> h;

    if (std::cin.fail())
    {
        std::cerr << "Error: could not read t or q from stdin\n";
        return 1;
    }

    // Perturb parameter
    std::cout.precision(17);
    // q = log10(q);
    // std::normal_distribution<double> distribution(q, stdev);
    // double q_perturbed = distribution(generator);
    std::lognormal_distribution<double> s_distribution(log(s), stdev1);
    std::normal_distribution<double> h_distribution(h, stdev2);
    double s_perturbed = s_distribution(generator);
    double h_perturbed = h_distribution(generator);
    // Print perturbed parameter
    std::cout << s_perturbed << " " << h_perturbed << std::endl;

    return 0;
}
