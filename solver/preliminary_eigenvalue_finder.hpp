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
	
    Eigenvalue locate_eigenvalue( double E_min, double E_max, const int nodes1, const int nodes2, const double a, const double b, const double h, const double eps );
    std::vector<Eigenvalue> findEigenvalues( const double E_min, const double E_max, const double a, const double b, const double h, double E_step = 1.0e-4, const int reduce_step = 20, const double eps = 1.0e-8 );

    std::vector<Eigenvalue> separate_eigenvalues( double E_min, double E_max, const int nodes_min, const int nodes_max, const double a, const double b, const double h, const double eps );

    void show_vector( std::string const & name, std::vector<Eigenvalue> const & v );

private:
	Equations * equations;
	int NPoints;
};
