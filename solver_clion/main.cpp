#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <map>
#include <gsl/gsl_sf_legendre.h>

#include "constants.hpp"
#include "equations.hpp"
#include "preliminary_eigenvalue_finder.hpp"
#include "precise_eigenvalue_finder.hpp"

#define ENABLE_WAVEFUNCTION_CALCULATION
//#undef ENABLE_WAVEFUNCTION_CALCULATION

const int NPoints = 8000;

void fixNodeCount( std::vector<Eigenvalue> & eigs );

std::vector<Eigenvalue> calculate_eigenvalues( Equations & equations,
        std::map<double, std::pair<double, double>> const & energy_dict,
        int NPoints, int channels, double E_min, double E_max,
        int i_match_intervals, double eps );

std::vector<Eigen::VectorXd> calculate_eigenfunction( Equations & equations,
        double E, int i_match, int NPoints, int channels, double & a, double & b, double & h );

int main( int argc, char * argv[] )
{
    std::vector<double> shift(3);

    int J, M, channels;
    if (argc < 7) {
        std::cerr << "Program expects (int) J, (int) M, (int) channels" << std::endl;
        exit(1);
    } else {
        J = atoi(argv[1]);
        M = atoi(argv[2]);
        channels = atoi(argv[3]);
        shift[0] = std::stod( argv[4] );
        shift[1] = std::stod( argv[5] );
        shift[2] = std::stod( argv[6] );

        std::cout << "(main) accepted arguments. J: " << J << ", M: " << M << "; channels: " << channels
                  << "; NPoints: " << NPoints << "; shifts: [" << shift[0] << ", " << shift[1] << ", " << shift[2] << "]" << std::endl;
    }

    std::cout << std::fixed << std::setprecision(16);

    // energy interval to search eigenvalues in
    const double E_min = -150.0 / constants::HTOCM;
    const double E_max = -0.1 / constants::HTOCM;
    const double eps = 1.0e-2; // relative precision of preliminary eigenvalue

    // interval to search turning points in  
    const double x_lb = 5.0;
    const double x_rb = 32.0;
    const double eps_tp = 1.0e-3; // absolute precision of turning point

    int energy_intervals = 10; // number of energy intervals [lg(-E_min), lg(-E_max)] to be divided into
    int i_match_intervals = 20; // number of dots condisered to be i_match

    std::vector<Eigenvalue> tmp, preigs;

    Equations equations(channels, NPoints, shift);
    equations.setAngularMomentum(J, M);

    std::map<double, std::pair<double, double>> energy_dict = equations.create_energy_dict(E_min, E_max,
                                                                                           energy_intervals, x_lb, x_rb,
                                                                                           eps_tp);

    std::cout << "Turning points corresponding to different energies: " << std::endl;
    bool all_tp_zero = true; // suppose turning point cannot be found
    for (auto it = energy_dict.begin(); it != energy_dict.end(); ++it) 
    {
        std::cout << it->first*constants::HTOCM << " (cm-1); (" << it->second.first << ", " << it->second.second << ")" << std::endl;
        // check if there exists at least one pair of non-zero turning points
        if ((it->second.first != 0.0) && (it->second.second != 0.0))
            all_tp_zero = false;
    }
    if (all_tp_zero)
    {
        std::cerr << "All turning points are zero! There are no bound levels. Check channels." << std::endl;
        exit( 1 );
    }

    tmp = calculate_eigenvalues( equations, energy_dict, NPoints, channels, E_min, E_max, i_match_intervals, eps );
    if ( !tmp.empty() && tmp[0].get_node_count_min() != 0 )
    {
        std::cerr << "Node count min for first eigenvalue is not zero! Exiting..." << std::endl;
        exit( 1 );
    }

    preigs.insert( preigs.end(), tmp.begin(), tmp.end() );

    std::cout << "List of sorted eigenvalues: " << std::endl;
    for ( size_t k = 0; k < preigs.size(); ++k ) {
        std::cout << k << " " << preigs[k].get_value()*constants::HTOCM << std::endl;
    }

	return 0;
}

std::vector<Eigenvalue> calculate_eigenvalues( Equations & equations, std::map<double, std::pair<double, double>> const & energy_dict,
        const int NPoints, const int channels, const double E_min, const double E_max, const int i_match_intervals, const double eps  )
{
    std::vector<Eigenvalue> preigs;
    PreliminaryEigenvalueFinder pef( &equations, energy_dict );
    std::vector<Eigenvalue> eigs = pef.findEigenvalues( E_min, E_max );

    std::cout << "Preliminary eigenvalues: " << std::endl;
    for ( size_t k = 0; k < eigs.size(); ++k )
    {
        eigs[k] = pef.convergeToEigenvalue( eigs[k], std::abs(eigs[k].get_value()) * eps );
        std::cout << k << " " << eigs[k].get_value()*constants::HTOCM << " (in interval: " << 
            eigs[k].get_min()*constants::HTOCM << ", " << eigs[k].get_max()*constants::HTOCM << ")" << std::endl;
    }

    PreciseEigenvalueFinder pref( &equations, energy_dict );
    pref.reserve_space_for_search( i_match_intervals );

    int i_match;
    std::cout << "Precise eigenvalues: " << std::endl;
    for ( size_t k = 0; k < eigs.size(); ++k )
    {
        i_match = pref.effecient_i_match( eigs[k].get_min(), eigs[k].get_max(), NPoints, i_match_intervals );
        std::cout << "(main) i_match: " << i_match << std::endl;
        if ( i_match < 0 )
            continue;

        double a, b, h;
        double eig = pref.precise_eigenvalue_calculation( i_match, eigs[k].get_min(), eigs[k].get_max(), a, b, h );
           
        preigs.push_back( eigs[k] );
        preigs.back().set_value( eig );
        //std::cout << "preigs.back().get_value(): " << preigs.back().get_value() * constants::HTOCM << std::endl;

        #ifdef ENABLE_WAVEFUNCTION_CALCULATION
        std::cout << "Calculating wavefunction..." << std::endl;
        std::vector<Eigen::VectorXd> eifunc = calculate_eigenfunction( equations, preigs.back().get_value(), i_match,
                NPoints, channels, a, b, h );

        std::ofstream ofs( "eifuncs/eifunc" + std::to_string(k) + ".txt" );
        ofs << std::setprecision(16);
        ofs << "# Energy: " << preigs.back().get_value() * constants::HTOCM << std::endl;
        ofs << "# J: " << equations.get_J() << std::endl;
        ofs << "# M: " << equations.get_M() << std::endl;
        ofs << "# Lmin: " << equations.get_M_abs() << std::endl;
        ofs << "# Lmax: " << equations.get_M_abs() + channels << std::endl;
        ofs << "# grid size: " << NPoints << std::endl;
        ofs << "# h: " << h << std::endl;
        for ( size_t i = 0; i < eifunc.size(); ++i )
        {
            ofs << a + i * h << " ";
            for ( int j = 0; j < eifunc[k].size(); ++j )
            {
                ofs << eifunc[i](j) << " ";
            }
            ofs << std::endl;
        }
        #endif

        //std::cout << preigs.back() << std::endl;
    }

    return preigs;
}

std::vector<Eigen::VectorXd> calculate_eigenfunction( Equations & equations, const double E, const int i_match, const int NPoints,
        const int channels, double & a, double & b, double & h )
{
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "(calculate_eigenfunction) Eigenvalue: " << E * constants::HTOCM << std::endl;

    Eigen::MatrixXd Rm, Rmp1;
    //std::pair<double, double> tp = equations.interpolate( E, energy_dict );
    //equations.calculate_boundaries( tp, &a, &b, &h );

    std::cout << "(eigenfunction) a: " << a << "; b: " << b << "; i_match: " << i_match << std::endl;
    equations.propagateForwardToMatch( E, a, h, i_match, Rm, true );
    equations.propagateBackwardToMatch( E, b, h, i_match, Rmp1, true );
    Eigen::MatrixXd diff = Rm - Rmp1.inverse();

    //std::cout << "diff: " << std::endl << diff << std::endl << std::endl;
    std::cout << "diff determinant: " << diff.determinant() << std::endl;

    Eigen::FullPivLU<Eigen::MatrixXd> lu( diff );
    // Be careful with this threshold value; may cause SEGFAULT; decrease if happens (?)
    lu = lu.setThreshold( 1.0e-6 );

    Eigen::MatrixXd ker = lu.kernel();
    std::cout << "kernel: " << std::endl << ker << std::endl;

    std::vector<Eigen::VectorXd> eigfunc(NPoints);
    for ( int k = 0; k < NPoints; ++k )
        eigfunc[k] = Eigen::VectorXd::Zero(channels);

    Eigen::VectorXd fn_prev, fn, fn_next, psin;
    Eigen::VectorXd fn_matching_point;

    fn_matching_point = ker;
    fn_next = fn_matching_point;
    fn_prev = fn_matching_point;

    //std::cout << "---------------------------------------------" << std::endl << std::endl;
    //std::cout << "(calculate_eigenfunction)" << std::endl;

    for ( int k = i_match - 1;  k > 0; --k, fn_next = fn )
    {
        //std::cout << "k: " << k << "; Rinv: " << std::endl << equations.Rinv_vector[k] << std::endl;
        fn = equations.Rinv_vector[k] * fn_next;
        psin = equations.Winv_vector[k] * fn;
        eigfunc[k] = psin;

        //std::cout << "eigfunc[" << k << "]: " << std::endl << eigfunc[k] << std::endl;
    }

    for ( int k = i_match + 1; k < NPoints - 1; ++k, fn_prev = fn )
    {
        //std::cout << "k: " << k << "; Rinv: " << std::endl << Rinv_vector[k] << std::endl;
        fn = equations.Rinv_vector[k] * fn_prev;
        psin = equations.Winv_vector[k] * fn;
        eigfunc[k] = psin;

        //std::cout << "eigfunc[" << k << "]: " << std::endl << eigfunc[k] << std::endl;
    }

    psin = equations.Winv_vector[i_match] * fn_matching_point;
    eigfunc[i_match] = psin;

    //std::cout << "eigfunc[" << i_match << "]: " << std::endl << eigfunc[i_match] << std::endl;

    return eigfunc;
}
