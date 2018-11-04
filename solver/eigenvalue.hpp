#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "./parity.hpp"

struct Eigenvalue
{
	Eigenvalue();
    
    Eigenvalue( Parity parity, double value, double min, double max, int node_count_min, int node_count_max  ); 
	Eigenvalue( Parity parity, double min, double max, int node_count_min, int node_count_max ); 
    Eigenvalue( Parity parity, double value, int node_count_min, int node_count_max );

    void setAngularMomentum( int J, int M );

    friend std::ostream& operator<< (std::ostream &o, const Eigenvalue & eig);
    friend std::ifstream& operator>>( std::ifstream & ifs, Eigenvalue & eig );

    bool operator<( const Eigenvalue & other ) const
    {
        return (this->J < other.get_J()) || 
              ((this->J == other.get_J()) && (this->M < other.get_M())) ||
              ((this->J == other.get_J()) && (this->M == other.get_M()) && (this->node_count_min < other.get_node_count_min())) ||
              ((this->J == other.get_J()) && (this->M == other.get_M()) && (this->node_count_min == other.get_node_count_min()) && (this->value < other.get_value())); 
    }

    //bool operator<( const Eigenvalue & other ) const
    //{
        //return value < other.value;
    //}

    void set_J( int J ) { this->J = J; }
    void set_M( int M ) { this->M = M; }
    void set_min( double min ) { this->min = min; }
    void set_max( double max ) { this->max = max; }
    void set_node_count_min( int node_count_min ) { this->node_count_min = node_count_min; }
    void set_node_count_max( int node_count_max ) { this->node_count_max = node_count_max; }
    void set_value( double value ) { this->value = value; } 
    void set_parity( Parity parity ) { this->parity = parity; }

    void update_value( ) { value = 0.5 * (min + max); }
    void dump_node_count() { node_count_min = 0; node_count_max = 0; } 

    int get_J() const { return J; }
    int get_M() const { return M; }
    int get_node_count_min() const { return node_count_min; }
    int get_node_count_max() const { return node_count_max; }
    
    double get_value() const { return value; }
    double get_min() const { return min; }
    double get_max() const { return max; }

private:
    Parity parity;
    double value;

    int J, M;
    double min, max;
    int node_count_min, node_count_max;

    int degeneracy;

    std::string type;
};
