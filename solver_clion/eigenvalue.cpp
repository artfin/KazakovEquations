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

Eigenvalue::Eigenvalue(  double min, double max, int node_count_min, int node_count_max ) :
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

void Eigenvalue::setAngularMomentum( int J, int M )
{
    this->J = J;
    this->M = M;
}
