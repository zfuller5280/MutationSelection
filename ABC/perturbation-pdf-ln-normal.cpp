#include <string>
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
double normal_pdf(long double mu, double sigma, long double x)
{
    // mu = log10(mu);
    // // return log10(exp(1)) * exp( - (log10(x) - mu) * (log10(x) - mu) / ( 2.0 * sigma * sigma ) ) /
    // //              (x * sigma * sqrt( 2.0 * M_PI ));
    // return (log10(exp(1)) /  (x * sigma * sqrt(2.0 * M_PI))) * exp(-0.5 * pow((log10(x)-mu)/sigma, 2));
    using boost::math::lognormal_distribution;
    lognormal_distribution<>my_dist(log(mu), sigma);
    long double prob = pdf(my_dist, (x));
    return prob;
}

int main(int argc, char *argv[])
{
    // Process arguments
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " STDEV\n"
            "Return probability densities of normal distributions centered \n"
            "on given parameters for the perturbed parameter \n"
            "and with standard deviation STDEV\n";
        return 1;
    }

    double stdev = std::stod(argv[1]);

    // Read generation t and parameter from stdin
    int t;
    std::cin >> t;

    long double q_perturbed;
    std::cin >> q_perturbed;

    if (std::cin.fail())
    {
        std::cerr << "Error: could not read t or q from stdin\n";
        return 1;
    }

    // Return normal_pdf for every parameter
    std::cout.precision(10);
    std::vector<double> q_vec;
    long double q;
    for (std::cin >> q; !std::cin.eof(); std::cin >> q)
    {
        if (std::cin.fail())
        {
            std::cerr << "Error: could not read parameter from stdin\n";
            return 1;
        }

        q_vec.push_back(q);
    }

    for (const double& q : q_vec)
    {
        long double pdf = normal_pdf(q, stdev, q_perturbed);
        if (pdf <= 1e-20){
          pdf = 0.;
        }
        std::cout << pdf << std::endl;
    }

    return 0;
}
