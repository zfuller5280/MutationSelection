#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <boost/math/distributions/beta.hpp>

/** Return probabilility density function of beta distribution.
 *
 * @param mu  mean of normal distribution.
 * @param sigma  standard deviation of normal distribution.
 * @param x  position to evaluate distribution on.
 *
 * @return probability density of normal distribution at x.
 */


using namespace boost::math;

double beta_pdf(double expec, double scale, double x)
{
    double alpha = expec * scale;
    double beta = (1-expec) * scale;
    beta_distribution<> d(alpha, beta);
    return pdf(d, x);
}

int main(int argc, char *argv[])
{
    // Process arguments
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " SCALE\n"
            "Return probability densities of normal distributions centered \n"
            "on given parameters for the perturbed parameter \n"
            "and with standard deviation SCALE\n";
        return 1;
    }

    double scale = std::stod(argv[1]);

    // Read generation t and parameter from stdin
    int t;
    std::cin >> t;

    double q_perturbed;
    std::cin >> q_perturbed;

    if (std::cin.fail())
    {
        std::cerr << "Error: could not read t or q from stdin\n";
        return 1;
    }

    // Return normal_pdf for every parameter
    std::cout.precision(17);
    std::vector<double> q_vec;
    double q;
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
        double pdf = beta_pdf(q, scale, q_perturbed);
        std::cout << pdf << std::endl;
    }

    return 0;
}
