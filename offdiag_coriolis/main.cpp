#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <functional>
#include <chrono>

#include <gsl/gsl_sf_coupling.h>
#include "wavefunction.hpp"
#include "filelister.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// taken from Manolopoulos Thesis
const double MU = 34505.15; 
// Hartree to inverse cm
const double HTOCM = 219475.797;

double LAMBDAP( const int A, const int B )
{
    return std::sqrt( A*(A+1) - B*(B+1) );
}

double LAMBDAM( const int A, const int B )
{
    return std::sqrt( A*(A+1) - B*(B-1) );
}

double channel_coriolis_matrix_element( std::vector<double> const& Rgrid, Channel const& ch1, Channel const& ch2,
        const int J, const int M, const int Mprime )
{
    int L = ch1.get_L();
    int Lprime = ch2.get_L();

    if ( L != Lprime ) return 0.0;

    if ( M != Mprime - 1 && M != Mprime + 1 )
        return 0.0;

    Channel product = ch1 * ch2;
    double radial_integral = integrate_simpson_r2(Rgrid, product);

    if ( M == Mprime + 1 ) {
        return 1.0 / (2.0 * MU) * radial_integral * LAMBDAP(J, Mprime) * LAMBDAP(Lprime, Mprime);
    }
    if ( M == Mprime - 1 ) {
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

double wf_coriolis_matrix_element( Wavefunction const& wf1, Wavefunction const& wf2 )
{
    double result = 0.0;
    
    int J = wf1.get_J();
    int Jprime = wf1.get_J();
    if ( J != Jprime ) return 0.0;

    int M = wf1.get_M();
    int Mprime = wf2.get_M();

    std::vector<double> const& Rgrid = wf1.get_grid();

    for ( Channel const& ch1 : wf1.get_channels() )
        for ( Channel const& ch2 : wf2.get_channels() )
            result += channel_coriolis_matrix_element(Rgrid, ch1, ch2, J, M, Mprime);

    return result * HTOCM;
}

double first_order(std::vector<Wavefunction> const& wavefunctions, const int n)
{
    return wf_coriolis_matrix_element(wavefunctions[n], wavefunctions[n]);
}

double second_order( std::vector<Wavefunction> const& wavefunctions, const int n)
// second-order perturbation correction for the n-th level by the Coriolis interaction
{
    size_t size_ = wavefunctions.size();
    double matrix_element, result = 0.0;
    double energy_n = wavefunctions[n].get_energy();

    for ( size_t k = 0; k < size_; ++k )
    {
        if ( k == n ) continue;
        matrix_element = wf_coriolis_matrix_element(wavefunctions[k], wavefunctions[n]); 
        result += matrix_element * matrix_element / (energy_n - wavefunctions[k].get_energy());
    }

    return result;
}

double third_order( std::vector<Wavefunction> const& wavefunctions, const int n )
// third-order perturbation correction to the energy of the n-th level by Coriolis interaction
// in the special case when first-order correction V_nn = 0
{
    size_t size_ = wavefunctions.size();
    double nk3, k3k2, k2n;
    double E_nk3, E_nk2;
    double E_n = wavefunctions[n].get_energy();

    double result = 0.0;
    for ( size_t k3 = 0; k3 < size_; ++k3 )
    {
        if ( k3 == n ) continue;
        nk3 = wf_coriolis_matrix_element( wavefunctions[n], wavefunctions[k3] );
        E_nk3 = E_n - wavefunctions[k3].get_energy();

        for ( size_t k2 = 0; k2 < size_; ++k2 )
        {
            if ( k2 == n ) continue;
            k3k2 = wf_coriolis_matrix_element(wavefunctions[k3], wavefunctions[k2]);
            k2n = wf_coriolis_matrix_element(wavefunctions[k2], wavefunctions[n]);
            E_nk2 = E_n - wavefunctions[k2].get_energy();

            result += nk3 * k3k2 * k2n / E_nk3 / E_nk2;
        }
    }

    return result;
}

double fourth_order( std::vector<Wavefunction> const& wavefunctions, const int n)
// fourth-order perturbation correction to the energy of the n-th level by Coriolis interaction
// in the special case when first-order correction V_nn = 0
{
    size_t size_ = wavefunctions.size();
    double nk4, k4k3, k3k2, nk2;
    double E_nk4, E_nk3, E_nk2;
    double E_n = wavefunctions[n].get_energy();

    double term1 = 0.0;
    for ( size_t k4 = 0; k4 < size_; ++k4 )
    {
        if ( k4 == n ) continue;
        nk4 = wf_coriolis_matrix_element( wavefunctions[n], wavefunctions[k4] );
        E_nk4 = E_n - wavefunctions[k4].get_energy();

        for ( size_t k3 = 0; k3 < size_; ++k3 )
        {
            if ( k3 == n ) continue;
            k4k3 = wf_coriolis_matrix_element( wavefunctions[k4], wavefunctions[k3] );
            E_nk3 = E_n - wavefunctions[k3].get_energy();

            for ( size_t k2 = 0; k2 < size_; ++k2 )
            {
                if ( k2 == n ) continue;
                k3k2 = wf_coriolis_matrix_element( wavefunctions[k3], wavefunctions[k2] );
                nk2 = wf_coriolis_matrix_element( wavefunctions[n], wavefunctions[k2] );
                E_nk2 = E_n - wavefunctions[k2].get_energy();
            
                term1 += nk4 * k4k3 * k3k2 * nk2 / E_nk4 / E_nk3 / E_nk2;
            }
        }

        //if ( k4 % 5 == 0 ) std::cout << "k4: " << k4 << std::endl;
    }

    double term2 = 0.0;
    for ( size_t k4 = 0; k4 < size_; ++k4 )
    {
        if ( k4 == n ) continue;
        nk4 = wf_coriolis_matrix_element( wavefunctions[n], wavefunctions[k4] );
        E_nk4 = E_n - wavefunctions[k4].get_energy();
        term2 += nk4 * nk4 / E_nk4 / E_nk4; 
    }

    double p2 = second_order( wavefunctions, n );
    return term1 - p2 * term2;
}

void print_results( std::vector<double> const& energies, const int n, std::vector<double> const& corrections, std::vector<double> const& eps_vector )
{
    std::cout << std::fixed << std::setprecision(9);
    double E0 = energies[n];
    double E;

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

double timer( std::function<double(std::vector<Wavefunction> const&, const int)> const& f, const std::vector<Wavefunction> & wavefunctions,
        int n, double & result )
{
    auto start = std::chrono::high_resolution_clock::now();
    result = f(wavefunctions, n);
    auto end = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
}

double calculate_m2_element( std::vector<Wavefunction> const& wavefunctions, const Wavefunction & wf1, const Wavefunction & wf2 )
// wavefunctions -- wavefunctions spanning non-degenerate subspace
// wf1, wf2 -- functions with
{
    double element = 0.0;
    double n1p, pn2;
    double E_n = wf1.get_energy(); // energy of wf2 is the same

    size_t size_ = wavefunctions.size();
    for ( size_t p = 0; p < size_; ++p )
    {
        n1p = wf_coriolis_matrix_element(wf1, wavefunctions[p]);
        pn2 = wf_coriolis_matrix_element(wavefunctions[p], wf2);
        element += n1p * pn2 / (E_n - wavefunctions[p].get_energy());
        //std::cout << "p: " << p << "; En - Ep: " << E_n - wavefunctions[p].get_energy() << std::endl;
    }

    return element;
}

double degenerate_second_order( std::vector<Wavefunction> const& wavefunctions, const int n1, const int n2 )
{
    //std::cout << "For second-order degenerate perturbation theory selected a pair of levels " << n1 << ", " << n2 << std::endl;
    //std::cout << "E1: " << wavefunctions[n1].get_energy() << "; E2: " << wavefunctions[n2].get_energy() << std::endl;
    //std::cout << "M1: " << wavefunctions[n1].get_M() << "; M2: " << wavefunctions[n2].get_M() << std::endl;

    // second-order 'mass' matrix (?)
    Eigen::Matrix2d M2 = Eigen::Matrix2d::Zero(2, 2);

    // vector of wavefunctions spanning non-degenerate subspace
    std::vector<Wavefunction> nd_wavefunctions = wavefunctions;
    nd_wavefunctions.erase(nd_wavefunctions.begin() + n1);
    nd_wavefunctions.erase(nd_wavefunctions.begin() + n2 - 1); // index is moving when we erase n1-th element...

    M2(0, 0) = calculate_m2_element(nd_wavefunctions, wavefunctions[n1], wavefunctions[n1]);
    M2(0, 1) = calculate_m2_element(nd_wavefunctions, wavefunctions[n1], wavefunctions[n2]);
    M2(1, 0) = calculate_m2_element(nd_wavefunctions, wavefunctions[n2], wavefunctions[n1]);
    M2(1, 1) = calculate_m2_element(nd_wavefunctions, wavefunctions[n2], wavefunctions[n2]);

    //std::cout << "M2 matrix: " << std::endl << M2 << std::endl;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es;
    es.compute( M2 );
    Eigen::Vector2d eigs = es.eigenvalues();

    if ( eigs(0) != 0.0 ) return eigs(0);
    return eigs(1);

}

int main()
{
    std::cout << std::fixed << std::setprecision(15);

    std::vector<std::string> dirs = {"../eifuncs_j1m0_nch12/", "../eifuncs_j1m1_nch12/", "../eifuncs_j1m-1_nch12/"};
    std::vector<std::string> files;
    for ( std::string const& dir : dirs )
    {
        FileLister lister(dir);
        std::vector<std::string> const& tmp = lister.get_files();
        files.insert(files.end(), tmp.begin(), tmp.end());
    }

    std::vector<Wavefunction> wavefunctions;
    for ( std::string const& file : files )
    {
        Wavefunction wavefunction;
        wavefunction.read( file );
        wavefunction.normalize();
        wavefunctions.push_back(wavefunction);
    }

    auto count = [&](int n) {
        return std::count_if(wavefunctions.begin(), wavefunctions.end(), [=](const Wavefunction& wf) { return wf.get_M() == n; });
    };
    std::cout << "M = 0 wavefunctions: " << count(0) << std::endl;
    std::cout << "M = 1 wavefunctions: " << count(1) << std::endl;
    std::cout << "M = -1 wavefunctions: " << count(-1) << std::endl;
    std::cout << "Total number of wavefunctions: " << wavefunctions.size() << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl << std::endl;

    /*
     * Degenerate second order perturbation theory
    std::cout << std::setprecision(4);
    for ( int k = static_cast<int>(count(0)); k < static_cast<int>(count(0) + count(1)); ++k )
    {
        int n1 = k;
        int n2 = k + static_cast<int>(count(1));
        //std::cout << "n1: " << n1 << "; n2: " << n2 << std::endl;

        double p2 = degenerate_second_order(wavefunctions, n1, n2);
        double E_n = wavefunctions[n1].get_energy();

        std::cout << "E_n(1): " << E_n << "; E_n(2): " << E_n + p2 << std::endl;
    }
     */

   /*
    int n = 1; // level we are correcting
    double E0 = wavefunctions[n].get_energy();
    std::cout << "Correcting level: " << n << "; E: " << E0 << "; J: " << wavefunctions[n].get_J() <<
                 "; M: " << wavefunctions[n].get_M() << std::endl;

    double p1, p2, p3, p4;
    double time1 = timer(first_order, wavefunctions, n, p1);
    std::cout << "(1 order) time: " << time1 << " ms." << std::endl;

    double time2 = timer(second_order, wavefunctions, n, p2);
    std::cout << "(2 order) time: " << time2 << " ms." << std::endl;

    double time3 = timer(third_order, wavefunctions, n, p3);
    std::cout << "(3 order) time: " << time3 << " ms." << std::endl;

    double time4 = timer(fourth_order, wavefunctions, n, p4);
    std::cout << "(4 order) time: " << time4 << " ms." << std::endl;

    std::cout << "p1: " << p1 << std::endl;
    std::cout << "p2: " << p2 << std::endl;
    //std::cout << "p3: " << p3 << std::endl;
    //std::cout << "p4: " << p4 << std::endl;
    std::cout << "En(0) = " << E0 << std::endl;
    std::cout << "En(2) = " << E0 + p2 << std::endl;
    //std::cout << "En(4) = " << E0 + p2 + p4 << std::endl;
    */

    std::cout << std::setprecision(4);
    auto nlvl = count(0);
    for ( int n = 0; n < nlvl; ++n )
    {
        double E0 = wavefunctions[n].get_energy();
        double p2 = second_order(wavefunctions, n);
        double p4 = fourth_order(wavefunctions, n);
        std::cout << n << " " << E0 << " " << E0 + p2 << " " << E0 + p2 + p4 << std::endl;
    }

    return 0;
}
