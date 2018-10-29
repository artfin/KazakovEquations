#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <functional>

#include <Eigen/Dense>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "./equations.hpp"

class PreciseEigenvalueFinder
{
public:
    PreciseEigenvalueFinder( Equations * equations, int NPoints ) : equations(equations), NPoints(NPoints) { }
        
    double brent( std::function<double(double)> f, double xb1, double xb2 );
    double D( const double a, const double b, const double h, const int i_match, const double E );
    double precise_eigenvalue_calculation( const double a, const double b, const double h, const int i_match, const double E1, const double E2 );
    
    int find_i_match( const double a, const double b, const double h, const double E1, const double E2, const int parts = 10 );

private:
    Equations * equations;
    int NPoints;    
};
