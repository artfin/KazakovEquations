//
// Created by artfin on 25.02.19.
//

#include "wavefunction.hpp"

std::vector<std::string> Wavefunction::split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

std::string Wavefunction::strip_spaces( std::string str )
{
    str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
    return str;
}

void Wavefunction::parameter_parser( const std::vector<std::string>& parameter_lines )
{
    for ( std::string const& line : parameter_lines )
    {
        std::vector<std::string> words = split(line, ':');
        std::string token = strip_spaces(split(words[0], '#')[1]);
        double value = std::stod( words[1] );

        if ( token == "Energy" ) { Energy = value; }
        else if ( token == "J" ) { J = static_cast<int>(value); }
        else if ( token == "M" ) { M = static_cast<int>(value); }
        else if ( token == "Lmin" ) { Lmin = static_cast<int>(value); }
        else if ( token == "Lmax" ) { Lmax = static_cast<int>(value); }
        else if ( token == "gridsize" ) { grid_size = static_cast<size_t>(value); }
        else if ( token == "h" ) { h = value; }
        else {
            std::cerr << "Unknown token = " << token << "! Exiting..." << std::endl;
            exit( 1 );
        }
    }
}

void Wavefunction::allocate_channels( )
{
    for ( int L = Lmin; L < Lmax; ++L )
    {
        Channel ch(L, h, grid_size);
        channels.push_back(ch);
    }
    Rgrid.reserve( grid_size );
}

void Wavefunction::read( std::string const& filename )
{
    std::ifstream ifs( filename );
    std::string line;
    std::vector<double> v;
    double d, R;

    bool PARAMETERS_SET = false;
    std::vector<std::string> parameter_lines;

    while ( std::getline(ifs, line) )
    {
        std::istringstream iss(line);

        size_t found = line.find('#');
        if ( found != std::string::npos ) {
            parameter_lines.push_back( line );
        } else if ( !PARAMETERS_SET ){
            parameter_parser( parameter_lines );
            PARAMETERS_SET = true;
            allocate_channels();
            goto here;
        } else {
            here:
            iss >> R;
            while (iss >> d)
                v.push_back(d);

            Rgrid.push_back(R);
            push_back_channels(v);

            v.clear();
        }
    }
}

void Wavefunction::push_back_channels(std::vector<double> v)
{
    for (size_t k = 0; k < v.size(); ++k)
        channels[k].push_back(v[k]);
}

void Wavefunction::normalize()
{
    double term, norm_constant = 0.0;

    for ( Channel const& channel : channels )
    {
        Channel squared = channel * channel;
        term = integrate_simpson(Rgrid, squared);
        //std::cout << "term: " << term << std::endl;
        norm_constant += term;
    }
    norm_constant = std::sqrt(norm_constant);
    std::cout << "Normalization constant: " << norm_constant << std::endl;

    for ( Channel & channel : channels )
        channel /= norm_constant;
}

