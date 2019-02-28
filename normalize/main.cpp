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

double channel_matrix_element( Channel const& ch1, Channel const& ch2, const int lambda, const int M, const double h )
// gsl_sf_coupling_3j accepts the values of spins and their projections multiplied by two.
{
    Channel product = ch1 * ch2;
    double radial_integral = integrate_simpson(h, product);

    int L = ch1.get_L();
    int Lprime = ch2.get_L();

    return std::pow(-1.0, M) * std::sqrt( (2.0*L + 1.0) * (2.0*Lprime + 1.0) ) * \
           gsl_sf_coupling_3j(2*Lprime, 2*lambda,2*L, 0, 0, 0) * gsl_sf_coupling_3j(2*Lprime, 2*lambda, 2*L, -2*M, 0, 2*M) * \
           radial_integral;
}

double wf_matrix_element( Wavefunction const& wf1, Wavefunction const& wf2, const int lambda, const int M )
{
    double result = 0.0;
    double h = wf1.get_h();

    for ( Channel const& ch1 : wf1.get_channels() )
        for ( Channel const& ch2 : wf2.get_channels() )
            result += channel_matrix_element(ch1, ch2, lambda, M, h);

    return result;
}

std::vector<std::string> create_file_list( std::string const& folder, const size_t n )
{
    std::vector<std::string> file_list;
    for ( size_t k = 0; k < n; ++k )
        file_list.push_back( folder + "eifunc" + std::to_string(k) + ".txt" );
    return file_list;
}

double first_order( std::vector<Wavefunction> const& wavefunctions, const int n, const int lambda, const int M )
// first order perturbation correction for the n-th energy level by the lambda-perturbation
{
    return wf_matrix_element(wavefunctions[n], wavefunctions[n], lambda, M);
}

double second_order( std::vector<Wavefunction> const& wavefunctions, std::vector<double> const& energies, const int n, const int lambda, const int M)
// second order perturbation correction for the n-th energy level by the lambda-perturbation
{
    size_t size_ = wavefunctions.size();
    double matrix_element = 0.0, result = 0.0;

    for ( size_t k = 0; k < size_; ++k )
    {
        if ( k == n ) continue;
        matrix_element = wf_matrix_element(wavefunctions[k], wavefunctions[n], lambda, M);
        //std::cout << "matrix element: " << matrix_element << std::endl;
        result += matrix_element*matrix_element / (energies[n] - energies[k]);
    }

    return result;
}

double third_order( std::vector<Wavefunction> const& wavefunctions, std::vector<double> const& energies, const int n, const int lambda, const int M)
// third order perturbation correction for the n-th energy level by the lambda-perturbation
{
    size_t size_ = wavefunctions.size();
    double result = 0.0;
    double kn = 0.0, km = 0.0, mn = 0.0;
    double nn = wf_matrix_element(wavefunctions[n], wavefunctions[n], lambda, M);

    for ( size_t k = 0; k < size_; ++k )
    {
        if ( k == n ) continue;

        kn = wf_matrix_element(wavefunctions[k], wavefunctions[n], lambda, M);
        for ( size_t m = 0; m < size_; ++m )
        {
            if (m == n ) continue;
            mn = wf_matrix_element(wavefunctions[m], wavefunctions[n], lambda, M);
            km = wf_matrix_element(wavefunctions[k], wavefunctions[m], lambda, M);

            result += mn * km * kn / (energies[m] - energies[n]) / (energies[k] - energies[n]);
        }

        result -= nn * kn * kn / (energies[k] - energies[n]) / (energies[k] - energies[n]);
    }

    return result;
}

double fourth_order( std::vector<Wavefunction> const& wavefunctions, std::vector<double> const& energies, const int n, const int lambda, const int M)
// the fourth-order perturbation correction for the n-th energy level by the lambda-perturbation
{
    size_t size_ = wavefunctions.size();
    double nk4 = 0.0, k4k3 = 0.0, k3k2 = 0.0, k2n = 0.0, k3n = 0.0;
    double E_nk4 = 0.0, E_nk3 = 0.0, E_nk2 = 0.0;

    double nn = wf_matrix_element(wavefunctions[n], wavefunctions[n], lambda, M);

    double term1 = 0.0;
    for ( size_t k4 = 0; k4 < size_; ++k4 )
    {
        if (k4 == n) continue;
        nk4 = wf_matrix_element(wavefunctions[n], wavefunctions[k4], lambda, M);
        E_nk4 = energies[n] - energies[k4];

        for ( size_t k3 = 0; k3 < size_; ++k3 )
        {
            if (k3 == n) continue;
            k4k3 = wf_matrix_element(wavefunctions[k4], wavefunctions[k3], lambda, M);
            E_nk3 = energies[n] - energies[k3];

            for ( size_t k2 = 0; k2 < size_; ++k2 )
            {
                if (k2 == n) continue;
                k3k2 = wf_matrix_element(wavefunctions[k3], wavefunctions[k2], lambda, M);
                k2n = wf_matrix_element(wavefunctions[k2], wavefunctions[n], lambda, M);
                E_nk2 = energies[n] - energies[k2];

                term1 += nk4 * k4k3 * k3k2 * k2n / (E_nk2 * E_nk3 * E_nk4);
            }
        }
    }

    double E_n2 = second_order(wavefunctions, energies, n, lambda, M);
    double term2 = 0.0;
    for ( size_t k4 = 0; k4 < size_; ++k4 )
    {
        if ( k4 == n ) continue;
        nk4 = wf_matrix_element(wavefunctions[n], wavefunctions[k4], lambda, M);
        E_nk4 = energies[n] - energies[k4];

        term2 += E_n2 * nk4 * nk4 / E_nk4 / E_nk4;
    }

    double term3 = 0.0;
    for ( size_t k4 = 0; k4 < size_; ++k4 )
    {
        if ( k4 == n ) continue;
        E_nk4 = energies[n] - energies[k4];
        nk4 = wf_matrix_element(wavefunctions[n], wavefunctions[k4], lambda, M);

        for ( size_t k3 = 0; k3 < size_; ++k3 )
        {
            if (k3 == n) continue;
            E_nk3 = energies[n] - energies[k3];
            k4k3 = wf_matrix_element(wavefunctions[k4], wavefunctions[k3], lambda, M);
            k3n = wf_matrix_element(wavefunctions[k3], wavefunctions[n], lambda, M);
            term3 += 2.0 * nn * nk4 * k4k3 * k3n / E_nk3 / E_nk3 / E_nk4;
        }
    }

    double term4 = 0.0;
    for ( size_t k4 = 0; k4 < size_; ++k4 )
    {
        if (k4 == n) continue;
        E_nk4 = energies[n] - energies[k4];
        nk4 = wf_matrix_element(wavefunctions[n], wavefunctions[k4], lambda, M);

        term4 += nn*nn*nk4*nk4 / std::pow(E_nk4, 3);
    }

    std::cout << "term1: "<< term1 << std::endl;
    std::cout << "term2: "<< term2 << std::endl;
    std::cout << "term3: "<< term3 << std::endl;
    std::cout << "term4: "<< term4 << std::endl;

    return term1 - term2 - term3 + term4;
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

    std::vector<double> energies = read_energies("../energies_j0m0.txt");
    std::string folder = "../eifuncs_j0m0/";
    std::vector<std::string> file_list = create_file_list(folder, energies.size());
    for ( std::string const& filename : file_list )
        std::cout << filename << std::endl;

    int M = 0;
    int lambda = 2;

    int Lmax = 8; size_t grid_size = 5000;

    std::vector<Wavefunction> wavefunctions;
    for ( std::string const& filename : file_list )
    {
        Wavefunction wf(Lmax, grid_size);
        wf.read( filename );
        wf.set_h();
        wf.normalize();
        wavefunctions.push_back(wf);
    }

    int n = 0; // the level which we are correcting using PT
    double p1 = first_order(wavefunctions, n, lambda, M);
    std::cout << "First order: " << p1 << std::endl;
    double p2 = second_order(wavefunctions, energies, n, lambda, M);
    std::cout << "Second order: " << p2 << std::endl;
    double p3 = third_order(wavefunctions, energies, n, lambda, M);
    std::cout << "Third order: " << p3 << std::endl;
    //double p4 = fourth_order(wavefunctions, energies, n, lambda, M);
    //std::cout << "Fourth order: " << p4 << std::endl;

    std::vector<double> corrections{ p1, p2, p3, p4 };
    std::vector<double> eps_vector = { 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1 };
    print_results(energies, n, corrections, eps_vector);

    return 0;
}