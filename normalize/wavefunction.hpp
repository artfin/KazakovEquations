#pragma once

#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

#include "channel.hpp"

class Wavefunction
{
public:
    Wavefunction( const int Lmax, const size_t grid_size ) :
        Lmax(Lmax), grid_size(grid_size)
    {
        Rgrid.reserve( grid_size );
        for ( int k = 0; k < Lmax; ++k )
        {
            channels.emplace_back(k, grid_size);
        }
    }

    void read( std::string const& filename );
    void push_back_channels( std::vector<double> v );
    void normalize();
    void set_h();
    double get_h() const { return h; }

    std::vector<double> const& get_grid() const { return Rgrid; }
    Channel const& get_channel( int k ) const { return channels[k]; }

    std::vector<Channel> const& get_channels() const { return channels; }

private:
    int Lmax;
    size_t grid_size;
    double h;

    std::vector<Channel> channels;
    std::vector<double> Rgrid;
};

