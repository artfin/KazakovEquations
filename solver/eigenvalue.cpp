#include "./eigenvalue.hpp"

Eigenvalue::Eigenvalue()
{
}

Eigenvalue::Eigenvalue( double value, double min, double max, int node_count_min, int node_count_max  ) : 
		value(value), min(min), max(max), node_count_min(node_count_min), node_count_max(node_count_max) 
{
	degeneracy = node_count_max - node_count_min;
    type = "preliminary";
}

Eigenvalue::Eigenvalue( double min, double max, int node_count_min, int node_count_max ) :
	min(min), max(max), node_count_min(node_count_min), node_count_max(node_count_max) 
{
	value = 0.5 * (min + max);
	degeneracy = node_count_max - node_count_min;
    type = "preliminary";
}

Eigenvalue::Eigenvalue( double value, int node_count_min, int node_count_max ) :
    value(value), node_count_min(node_count_min), node_count_max(node_count_max)
{
    degeneracy = node_count_max - node_count_min;
    type = "precise";
}

std::ostream& operator<< (std::ostream &o, const Eigenvalue & eig)
{
    if ( eig.type == "preliminary" )
    	o << "(Eigenvalue) min: " << eig.min << "; max: " << eig.max << "; value: " << eig.value << "; NC min: " << eig.node_count_min << "; NC max: " << eig.node_count_max << "; degeneracy: " << eig.degeneracy;
    else if ( eig.type == "precise" )
        o << "(Eigenvalue) value: " << eig.value << "; NC min: " << eig.node_count_min << "; NC max: " << eig.node_count_max << "; degeneracy: " << eig.degeneracy; 

  	return o;
}

