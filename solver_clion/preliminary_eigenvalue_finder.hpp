#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <string>
#include <map>

#include <Eigen/Dense>

#include "./eigenvalue.hpp"
#include "./equations.hpp"

class PreliminaryEigenvalueFinder
{
public:
	PreliminaryEigenvalueFinder( Equations * equations, std::map<double, std::pair<double, double>> const& energy_dict ) :
		equations(equations), energy_dict(energy_dict)
    {
    }
	
    std::vector<Eigenvalue> findEigenvalues( double E_min, double E_max );
    std::vector<Eigenvalue> recursive_find( double E_min, double E_max, int nodes_min, int nodes_max );

    void show_vector( std::string const & name, std::vector<Eigenvalue> const& v );
    Eigenvalue convergeToEigenvalue( Eigenvalue const & e, double eps );

private:
	Equations * equations;
    const std::map<double, std::pair<double, double>> & energy_dict;
};
