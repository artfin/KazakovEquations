//
// Created by artfin on 25.02.19.
//
#include "channel.hpp"

Channel operator*( Channel lhs, const Channel & rhs )
{
    return lhs *= rhs;
}

double integrate_trapezoid( const double h, Channel const& ch )
{
    double sum = 0.0;
    size_t size_ = ch.get_grid_size();
    std::vector<double> const& y = ch.get_data();

    for (size_t k = 1; k < size_ - 1; ++k)
    {
        sum += 2.0 * y[k];
    }
    sum = sum + y[0] + y.back();

    return sum * h / 2.0;
}

double integrate_simpson( std::vector<double> const& Rgrid, Channel const& ch )
{
    double sum1 = 0.0; double sum2 = 0.0;
    size_t size_ = ch.get_grid_size();
    std::vector<double> const& y = ch.get_data();

    for ( size_t k = 1; k < size_ - 1; ++k )
    {
        if ( k % 2 == 0 ) {
            sum1 += y[k]; // / std::pow(Rgrid[k], 2);
        }
        else {
            sum2 += y[k]; // / std::pow(Rgrid[k], 2);
        }
    }

    double res = y[0] + y.back() + 2.0*sum1 + 4.0*sum2;
    return res * ch.get_h() / 3.0;
}

