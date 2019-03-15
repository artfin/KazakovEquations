#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

#include "./eigenvalue.hpp"

class EigenvalueDatabase 
{
public:
    EigenvalueDatabase( std::string db_filename ) : db_filename(db_filename) 
    {
    }

    bool file_exists( std::string const & filename );
    
    void read( );
    
    bool check_unique( const Eigenvalue & eig );
    
    void show();

    void insert( std::vector<Eigenvalue> const & eigs );

    void write();

private:
    std::string db_filename; 

    std::vector<Eigenvalue> db;
};

