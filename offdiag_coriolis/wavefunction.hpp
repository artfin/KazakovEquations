#pragma once

#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "channel.hpp"

class Wavefunction
{
public:
    explicit Wavefunction( )
    {
    }

    std::vector<std::string> split(const std::string& s, char delimiter);
    std::string strip_spaces( std::string str );
    void parameter_parser( const std::vector<std::string>& parameter_lines );

    void allocate_channels();

    void read( std::string const& filename );
    void push_back_channels( std::vector<double> v );

    void normalize();

    std::vector<double> const& get_grid() const { return Rgrid; }
    double get_h() const { return h; }
    size_t get_grid_size() const { return grid_size; }
    int get_J() const { return J; }
    int get_M() const { return M; }
    double get_energy() const { return Energy; }

    std::vector<Channel> const& get_channels() const { return channels; }
    Channel const& get_channel( int k ) const { return channels[k]; }

private:
    double Energy;
    int J;
    int M;
    int Lmin;
    int Lmax;
    size_t grid_size;
    double h;

    std::vector<Channel> channels;
    std::vector<double> Rgrid;
};

