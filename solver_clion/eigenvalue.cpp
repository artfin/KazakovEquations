#include "./eigenvalue.hpp"

Eigenvalue::Eigenvalue()
{
}

Eigenvalue::Eigenvalue( Parity parity, double value, double min, double max, int node_count_min, int node_count_max  ) : 
		parity(parity), value(value), min(min), max(max), node_count_min(node_count_min), node_count_max(node_count_max) 
{
	degeneracy = node_count_max - node_count_min;
    type = "preliminary";
}

Eigenvalue::Eigenvalue( Parity parity, double min, double max, int node_count_min, int node_count_max ) :
	parity(parity), min(min), max(max), node_count_min(node_count_min), node_count_max(node_count_max) 
{
	value = 0.5 * (min + max);
	degeneracy = node_count_max - node_count_min;
    type = "preliminary";
}

Eigenvalue::Eigenvalue( Parity parity, double value, int node_count_min, int node_count_max ) :
    parity(parity), value(value), node_count_min(node_count_min), node_count_max(node_count_max)
{
    degeneracy = node_count_max - node_count_min;
    type = "precise";
}

void Eigenvalue::setAngularMomentum( int J, int M )
{
    this->J = J;
    this->M = M;
}

std::ostream& operator<<( std::ostream& ofs, const Eigenvalue & eig )
{
    ofs << eig.J << " " << eig.M << " " << eig.node_count_min << " " << eig.node_count_max << " " << " " << eig.value * constants::HTOCM << " " << eig.parity; 
    return ofs;
}

std::ifstream& operator>>( std::ifstream& ifs, Eigenvalue & eig )
{
    int J, M, node_count_min, node_count_max;
    Parity p;
    double value;

    ifs >> J;
    eig.set_J( J );

    ifs >> M;
    eig.set_M( M );

    ifs >> node_count_min;
    eig.set_node_count_min( node_count_min ); 
    
    ifs >> node_count_max;
    eig.set_node_count_max( node_count_max );

    ifs >> value;
    eig.set_value( value );
    
    ifs >> p;
    eig.set_parity( p );

    return ifs; 
}

