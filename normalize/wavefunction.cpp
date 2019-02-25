//
// Created by artfin on 25.02.19.
//

#include "wavefunction.hpp"

void Wavefunction::read( std::string const& filename )
{
    std::ifstream ifs( filename );
    std::string line;
    std::vector<double> v;
    double d, R;

    while ( std::getline(ifs, line) )
    {
        std::istringstream iss(line);

        iss >> R;
        while (iss >> d)
            v.push_back( d );

        Rgrid.push_back( R );
        push_back_channels( v );

        v.clear();
    }
}

void Wavefunction::push_back_channels(std::vector<double> v)
{
    for (size_t k = 0; k < v.size(); ++k)
        channels[k].push_back(v[k]);
}

void Wavefunction::set_h()
{
    h = Rgrid[1] - Rgrid[0];
    std::cout << "h (bohrs): " << h << std::endl;
}

void Wavefunction::normalize()
{
    double term, norm_constant = 0.0;

    for ( Channel const& channel : channels )
    {
        Channel squared = channel * channel;
        term = integrate_simpson(h, squared);
        //std::cout << "term: " << term << std::endl;
        norm_constant += term;
    }
    norm_constant = std::sqrt(norm_constant);
    std::cout << "Normalization constant: " << norm_constant << std::endl;

    for ( Channel & channel : channels )
        channel /= norm_constant;
}

