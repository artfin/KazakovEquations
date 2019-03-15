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

#include "./equations.hpp"

class PreciseEigenvalueFinder
{
public:
    PreciseEigenvalueFinder( Equations * equations, std::map<double, std::pair<double, double>> const& energy_dict ) :
        equations(equations), energy_dict(energy_dict) { }
        
    double brent( std::function<double(double)> f, double xb1, double xb2 );
    double D( double a, double b, double h, int i_match, double E );
    double precise_eigenvalue_calculation( int i_match, double E1, double E2, double & a, double & b, double & h );
    
    int find_i_match( double E1, double E2, int NPoints, int parts = 100 );

    int effecient_i_match( double E1, double E2, int NPoints, int parts );

    void reserve_space_for_search( int parts );

private:
    Equations * equations;
    const std::map<double, std::pair<double, double>> & energy_dict;

    std::vector<Eigen::MatrixXd> Rm_vector1;
    std::vector<Eigen::MatrixXd> Rmp1_vector1;
    
    std::vector<Eigen::MatrixXd> Rm_vector2;
    std::vector<Eigen::MatrixXd> Rmp1_vector2;

    std::vector<double> D1_vector;
    std::vector<double> D2_vector;
};
