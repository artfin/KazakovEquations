#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <string>
#include <map>

#include <Eigen/Dense>

#include "./parity.hpp"
#include "./eigenvalue.hpp"
#include "./equations.hpp"

class PreliminaryEigenvalueFinder
{
public:
	PreliminaryEigenvalueFinder( Equations * equations, std::map<double, std::pair<double, double>> const & energy_dict, Parity parity ) :
		equations(equations), energy_dict(energy_dict), parity(parity) 
    {
    }
	
    std::vector<Eigenvalue> findEigenvalues( const double E_min, const double E_max );
    std::vector<Eigenvalue> recursive_find( const double E_min, const double E_max, const int nodes_min, const int nodes_max );

    void show_vector( std::string const & name, std::vector<Eigenvalue> const & v );
    Eigenvalue convergeToEigenvalue( Eigenvalue const & e, const double eps );

private:
	Equations * equations;
    std::map<double, std::pair<double, double>> energy_dict;
    Parity parity;
};
