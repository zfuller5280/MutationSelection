#include <string>
#include <iostream>
#include <math.h>

int main(int argc, char *argv[])
{
    // Process arguments
    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " S_LOW S_HIGH H_LOW H_HIGH\n"
            "Return probability density of log10-uniform distribution with \n"
            "lower bound S_LOW and upper bound S_HIGH "
            "for given parameter\n"
            "and uniformly distributed gamma with "
            "lower bound H_LOW and upper H_HIGH,\n"
            "for given parameter\n";
        return 1;
    }

    double s_low = std::stod(argv[1]);
    double s_high = std::stod(argv[2]);
    double h_low = std::stod(argv[3]);
    double h_high = std::stod(argv[4]);

    if (s_high <= s_low)
    {
        std::cerr << "Error: S_LOW must be strictly less than S_HIGH\n";
        return 1;
    }

    if (h_high <= h_low)
    {
        std::cerr << "Error: H_LOW must be strictly less than H_HIGH\n";
        return 1;
    }

    // Read parameter from stdin
    double s, h;
    std::cin >> s >> h;

    if (std::cin.fail())
    {
        std::cerr << "Error: could not read q from stdin\n";
        return 1;
    }

    // Check if q is within bounds
    std::cout.precision(17);
    if ( (s_low <= s) && (s <= s_high) && (h_low <= h) && (h <= h_high) )
    {
        std::cout << 1.0 / (s * (log10(s_high) - log10(s_low))) / (h_high - h_low)  << std::endl;
    }
    else
    {
        std::cout << 0 << std::endl;
    }

    return 0;
}
