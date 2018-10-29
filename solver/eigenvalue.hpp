#pragma once

#include <iostream>
#include <string>

struct Eigenvalue
{
	Eigenvalue();

    Eigenvalue( double value, double min, double max, int node_count_min, int node_count_max  ); 

	Eigenvalue( double min, double max, int node_count_min, int node_count_max ); 

    Eigenvalue( double value, int node_count_min, int node_count_max );

    friend std::ostream& operator<< (std::ostream &o, const Eigenvalue & eig);

    double value;
    
    double min;
    double max;

    int node_count_min;
    int node_count_max;

    int degeneracy;

    std::string type;
};
