#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cassert>

#include <gsl/gsl_sf_coupling.h>
#include "wavefunction.hpp"

// taken from Manolopoulos Thesis
const double MU = 34505.15; 
// Hartree to inverse cm
const double HTOCM = 219475.797;

std::vector<double> read_energies( std::string const& filename )
{
    std::ifstream ifs( filename );
    std::string line;

    std::vector<double> energies;
    double energy;

    while ( std::getline(ifs, line) )
    {
        std::cout << line << std::endl;
        std::istringstream iss(line);
        iss >> energy;
        energies.push_back( energy );
    }

    return energies;
}

double LAMBDAP( const int A, const int B )
{
    return std::sqrt( A*(A+1) - B*(B+1) );
}

double LAMBDAM( const int A, const int B )
{
    return std::sqrt( A*(A+1) - B*(B-1) );
}

double channel_coriolis_matrix_element( std::vector<double> const& Rgrid, Channel const& ch1, Channel const& ch2, const int J, const int M, const int Mprime )
{
    int L = ch1.get_L();
    int Lprime = ch2.get_L();

    //std::cout << "M: " << M << "; Mprime: " << Mprime << std::endl;
    //std::cout << "L: " << L << "; Lprime: " << Lprime << std::endl;

    if ( L != Lprime ) return 0.0;

    if ( M == Mprime + 1 )
    {
        Channel product = ch1 * ch2;
        double radial_integral = integrate_simpson_r2(Rgrid, product);
        //std::cout << "radial_integral: " << radial_integral << std::endl;
        return 1.0 / (2.0 * MU) * radial_integral * LAMBDAP(J, Mprime) * LAMBDAP(Lprime, Mprime);
    }
    if ( M == Mprime - 1 )
    {
        Channel product = ch1 * ch2;
        double radial_integral = integrate_simpson_r2(Rgrid, product);
        //std::cout << "radial_integral: " << radial_integral << std::endl;
        return 1.0 / (2.0 * MU) * radial_integral * LAMBDAM(J, Mprime) * LAMBDAM(Lprime, Mprime);
    }

    return 0.0;
}

/*
double channel_coriolis_matrix_element( std::vector<double> const& Rgrid, Channel const& ch1, Channel const& ch2, const int J, const int M, const int Mprime )
{
    int L = ch1.get_L();
    int Lprime = ch2.get_L();

    if ( L != Lprime ) return 0.0;
    if ( M != Mprime ) return 0.0;
    
    Channel product = ch1 * ch2;
    return integrate_simpson_r2(Rgrid, product);
}
*/

std::vector<std::string> create_file_list( std::string const& folder, const size_t n )
{
    std::vector<std::string> file_list;
    for ( size_t k = 0; k < n; ++k )
        file_list.push_back( folder + "eifunc" + std::to_string(k) + ".txt" );
    return file_list;
}

double wf_coriolis_matrix_element( Wavefunction const& wf1, Wavefunction const& wf2 )
{
    double result = 0.0;
    
    int J = wf1.get_J();
    int M = wf1.get_M();
    int Mprime = wf2.get_M();
    std::vector<double> const& Rgrid = wf1.get_grid();

    for ( Channel const& ch1 : wf1.get_channels() )
        for ( Channel const& ch2 : wf2.get_channels() )
            result += channel_coriolis_matrix_element(Rgrid, ch1, ch2, J, M, Mprime);

    return result;
}

double first_order(std::vector<Wavefunction> const& wavefunctions, const int n)
{
    return wf_coriolis_matrix_element(wavefunctions[n], wavefunctions[n]) * HTOCM;
}

double second_order( std::vector<Wavefunction> const& wavefunctions, std::vector<double> const& energies, const int n)
// second order perturbation correction for the n-th level by the Coriolis interaction
{
    size_t size_ = wavefunctions.size();
    double matrix_element = 0.0, result = 0.0;

    for ( size_t k = 0; k < size_; ++k )
    {
        if ( k == n ) continue;
        matrix_element = wf_coriolis_matrix_element(wavefunctions[k], wavefunctions[n]); 
        //std::cout << "matrix element: " << matrix_element << std::endl;
        result += matrix_element*matrix_element / (energies[n] - energies[k]);
    }

    return result * std::pow(HTOCM, 2);
}

double fourth_order( std::vector<Wavefunction> const& wavefunctions, std::vector<double> const& energies, const int n)
{
    size_t size_ = wavefunctions.size();
    double matrix_element= 0.0, result = 0.0;
    double nk4 = 0.0, k4k3 = 0.0, k3k2 = 0.0, nk2 = 0.0;
    double E_nk4 = 0.0, E_nk3 = 0.0, E_nk2 = 0.0;

    double term1 = 0.0;
    for ( size_t k4 = 0; k4 < size_; ++k4 )
    {
        if ( k4 == n ) continue;
        nk4 = wf_coriolis_matrix_element( wavefunctions[n], wavefunctions[k4] );
        E_nk4 = energies[n] - energies[k4];

        for ( size_t k3 = 0; k3 < size_; ++k3 )
        {
            if ( k3 == n ) continue;
            k4k3 = wf_coriolis_matrix_element( wavefunctions[k4], wavefunctions[k3] );
            E_nk3 = energies[n] - energies[k3];

            for ( size_t k2 = 0; k2 < size_; ++k2 )
            {
                if ( k2 == n ) continue;
                k3k2 = wf_coriolis_matrix_element( wavefunctions[k3], wavefunctions[k2] );
                nk2 = wf_coriolis_matrix_element( wavefunctions[n], wavefunctions[k2] );
                E_nk2 = energies[n] - energies[k2];
            
                term1 += nk4 * k4k3 * k3k2 * nk2 / E_nk4 / E_nk3 / E_nk2;
            }
        }
        std::cout << "k4: " << k4 << std::endl;
    }

    double term2 = 0.0;
    for ( size_t k4 = 0; k4 < size_; ++k4 )
    {
        if ( k4 == n ) continue;
        nk4 = wf_coriolis_matrix_element( wavefunctions[n], wavefunctions[k4] );
        E_nk4 = energies[n] - energies[k4];
        term2 += nk4 * nk4 / E_nk4 / E_nk4; 
    }

    double p2 = second_order( wavefunctions, energies, n );
    return term1 * std::pow(HTOCM, 4) - p2 * term2 * std::pow(HTOCM, 2);
}

void print_results( std::vector<double> const& energies, const int n, std::vector<double> const& corrections, std::vector<double> const& eps_vector )
{
    std::cout << std::fixed << std::setprecision(9);
    double E0 = energies[n];
    double E = E0;

    std::cout << "Epsilon \t corrections\n";
    std::cout << 0.0 << "\t" << E0 << std::endl;
    for ( double epsilon : eps_vector )
    {
        std::cout << epsilon << "\t";
        E = E0;

        for ( size_t k = 0; k < corrections.size(); ++k )
        {
            E += corrections[k] * std::pow(epsilon, k + 1);
            std::cout << E << "\t";
        }
        std::cout << "\n";
    }
}

int main()
{
    std::cout << std::fixed << std::setprecision(15);

    std::vector<double> energies_j1m0 = read_energies("./energies_j1m0.txt");
    std::vector<double> energies_j1m1 = read_energies("./energies_j1m1.txt");
    
    std::vector<double> energies = energies_j1m0;
    energies.insert( energies.end(), energies_j1m1.begin(), energies_j1m1.end() ); // M == 1
    energies.insert( energies.end(), energies_j1m1.begin(), energies_j1m1.end() ); // M == -1

    std::string folder1 = "./eifuncs_j1m0/";
    std::string folder2 = "./eifuncs_j1m1/";
    std::vector<std::string> file_list1 = create_file_list(folder1, energies_j1m0.size());
    std::vector<std::string> file_list2 = create_file_list(folder2, energies_j1m1.size());
    std::cout << file_list1.size() << " files for J = 1, M = 0 to read." << std::endl;
    std::cout << file_list2.size() << " files for J = 1, M = 1 to read." << std::endl;

    int J = 1;
    int M = 0;
    int lambda = 2;

    int Lmax = 8; size_t grid_size = 5000;

    std::vector<Wavefunction> wavefunctions;
    for ( std::string const& filename : file_list1 )
    {
        Wavefunction wf( J, M, Lmax, grid_size );
        wf.read( filename );
        wf.set_h();
        wf.normalize();
        wavefunctions.push_back(wf);
    }
    M = 1;
    for ( std::string const& filename : file_list2 )
    {
        Wavefunction wf( J, M, Lmax, grid_size );
        wf.read( filename );
        wf.set_h();
        wf.normalize();
        wavefunctions.push_back( wf );
    }
    
    M = -1;
    for ( std::string const& filename : file_list2 )
    {
        Wavefunction wf( J, M, Lmax, grid_size );
        wf.read( filename );
        wf.set_h();
        wf.normalize();
        wavefunctions.push_back( wf );
    }

    std::cout << wavefunctions.size() << " wavefunctions read." << std::endl;

    int n = 1; // level we are correcting
    double p2 = second_order(wavefunctions, energies, n);
    double p4 = fourth_order(wavefunctions, energies, n);

    std::cout << "p2: " << p2 << std::endl;
    std::cout << "p4: " << p4 << std::endl;
    std::cout << "En(0) = " << energies[n] << std::endl;
    std::cout << "En(2) = " << energies[n] + p2 << std::endl;
    std::cout << "En(4) = " << energies[n] + p2 + p4 << std::endl;

    return 0;
}
