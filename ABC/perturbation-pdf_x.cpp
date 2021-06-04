#include <string>
#include <utility>
#include <vector>
#include <iostream>

#include <math.h>
#include <boost/math/distributions/lognormal.hpp>

/** Return probabilility density function of normal distribution.
 *
 * @param mu  mean of normal distribution.
 * @param sigma  standard deviation of normal distribution.
 * @param x  position to evaluate distribution on.
 *
 * @return probability density of normal distribution at x.
 */
double lognormal_pdf(long double mu, double sigma, long double x)
{
    using boost::math::lognormal_distribution;
    lognormal_distribution<>my_dist(log(mu), sigma);
    long double prob = pdf(my_dist, (x));
    return prob;
}

// Define normal pdf
double normal_pdf(double mu, double sigma, double x)
{

    return exp( - (x - mu) * (x - mu) / ( 2.0 * sigma * sigma ) ) /
                 sqrt( 2.0 * M_PI * sigma * sigma );
}

int main(int argc, char *argv[])
{
    // Process arguments
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " STDEV1 STDEV2\n"
            "Return probability densities of joint log-normal distribution centered \n"
            "on given parameters for the perturbed parameter \n"
            "and with standard deviation STDEV1\n"
            "and normal distribution centered on given h"
            "and with standard deviation STDEV2";
        return 1;
    }

    double stdev1 = std::stod(argv[1]);
    double stdev2 = std::stod(argv[2]);

    // Read generation t and parameter from stdin
    int t;
    std::cin >> t;

    long double s_perturbed, h_perturbed;
    std::cin >> s_perturbed >> h_perturbed;

    if (std::cin.fail())
    {
        std::cerr << "Error: could not read t or q from stdin\n";
        return 1;
    }

    // Return normal_pdf for every parameter
    std::cout.precision(10);
    std::vector<std::pair<long double,double>> sh_vec;;
    long double s;
    double h;
    for (std::cin >> s >> h; !std::cin.eof(); std::cin >> s >> h)
    {
        if (std::cin.fail())
        {
            std::cerr << "Error: could not read parameter from stdin\n";
            return 1;
        }

        sh_vec.emplace_back(s, h);
    }

    for (const auto& sh : sh_vec)
    {
        long double s = sh.first;
        double h = sh.second;
        long double pdf1 = lognormal_pdf(s, stdev1, s_perturbed);
        if (pdf1 <= 1e-20){
          pdf1 = 0.;
        }
        double pdf2 = normal_pdf(h, stdev2, h_perturbed);
        double pdf = pdf1 * pdf2;
        std::cout << pdf << std::endl;
    }

    return 0;
}
