#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <map>
#include <functional>

#include <Eigen/Dense>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "./parity.hpp"
#include "./equations.hpp"

class PreciseEigenvalueFinder
{
public:
    PreciseEigenvalueFinder( Equations * equations, std::map<double, std::pair<double, double>> const & energy_dict, Parity parity ) : 
        equations(equations), energy_dict(energy_dict), parity(parity) { }
        
    double brent( std::function<double(double)> f, double xb1, double xb2 );
    double D( const double a, const double b, const double h, const int i_match, const double E );
    double precise_eigenvalue_calculation( const int i_match, const double E1, const double E2 );
    
    int find_i_match( const double E1, const double E2, const int NPoints, const int parts = 100 );

    int effecient_i_match( const double E1, const double E2, const int NPoints, const int parts );

    void reserve_space_for_search( const int parts );

private:
    Equations * equations;
    std::map<double, std::pair<double, double>> energy_dict;
    Parity parity;

    std::vector<Eigen::MatrixXd> Rm_vector1;
    std::vector<Eigen::MatrixXd> Rmp1_vector1;
    
    std::vector<Eigen::MatrixXd> Rm_vector2;
    std::vector<Eigen::MatrixXd> Rmp1_vector2;

    std::vector<double> D1_vector;
    std::vector<double> D2_vector;
};
