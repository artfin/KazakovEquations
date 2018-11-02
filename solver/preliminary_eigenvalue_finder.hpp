#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <string>

#include <Eigen/Dense>

#include "./eigenvalue.hpp"
#include "./equations.hpp"

class PreliminaryEigenvalueFinder
{
public:
	PreliminaryEigenvalueFinder( Equations * equations, int NPoints ) :
		equations(equations), NPoints(NPoints) 
    {
    }
	
    std::vector<Eigenvalue> findEigenvalues( const double E_min, const double E_max, const double a, const double b, const double h );
    std::vector<Eigenvalue> recursive_find( const double E_min, const double E_max, const int nodes_min, const int nodes_max, const double a, const double b, const double h );

    void show_vector( std::string const & name, std::vector<Eigenvalue> const & v );
    Eigenvalue convergeToEigenvalue( Eigenvalue const & e, const double a, const double b, const double h, const double eps );

private:
	Equations * equations;
	int NPoints;
};
