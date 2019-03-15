#pragma once

#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include "./equations.hpp"

class EigenfunctionFinder
{
public:
    EigenfunctionFinder( Equations * equations, int channels, int NPoints ) : equations(equations), channels(channels), NPoints(NPoints) 
    {
    } 

    void calculate_eigenfunction( std::vector<Eigen::VectorXd> & eigfunc, const double E, const double a, const double b, const double h, int i_match );

private:
    Equations * equations;
    int channels;
    int NPoints;
};
