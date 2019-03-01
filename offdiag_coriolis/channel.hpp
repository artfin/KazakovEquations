#pragma once

#include <iostream>
#include <cmath>
#include <vector>

class Channel
{
public:
    Channel( int L, double h, size_t grid_size ) : L(L), h(h), grid_size(grid_size)
    {
        data.reserve( grid_size );
    }

    int get_L() const { return L; }
    double get_h() const { return h; }
    size_t get_grid_size() const { return grid_size; }

    void push_back( double d )
    {
        data.push_back( d );
    }

    std::vector<double> & get_data() { return data; }
    std::vector<double> const& get_data() const { return data; }

    Channel& operator*= (const Channel& other)
    {
        std::vector<double> const& other_data = other.get_data();
        size_t size_ = other_data.size();

        for ( size_t k = 0; k < size_; ++k )
            data[k] *= other_data[k];

        return *this;
    }

    Channel& operator/=( const double d )
    {
        std::vector<double> & data = get_data();
        size_t size_ = data.size();

        for ( size_t k = 0; k < size_; ++k )
            data[k] /= d;

        return *this;
    }

private:
    int L;
    double h;
    size_t grid_size;

    std::vector<double> data;
};

Channel operator*( Channel lhs, const Channel & rhs );
double integrate_simpson( std::vector<double> const& Rgrid, Channel const& ch );
double integrate_simpson_r2( std::vector<double> const& Rgrid, Channel const& ch );


