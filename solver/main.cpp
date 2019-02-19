#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <map>

#include "./constants.hpp"

#include "./parity.hpp"
#include "./equations.hpp"
#include "./preliminary_eigenvalue_finder.hpp"
#include "./precise_eigenvalue_finder.hpp"
//#include "./eigenfunction_finder.hpp"

#include <gsl/gsl_sf_legendre.h>

#include "./eigenvalue_database.hpp"

// TO DO LIST:
// 0) discern parity of eigenvalues -- done
// 1) same automatic deduction of value of 'eps' in PreliminaryEigenvalueFinder::convergeToEigenvalue    ?!
// 
// 2) automatic update of interval bounds in PreliminaryEigenvalueFinder::findEigenvalues():  -- done
//      a = a(E), b = b(E)
// 2a) turning points depend on parity     
//
// 3) database to save eigenvalue into -- done 
// 3a) add J, M fields to eigenvalue structure -- done; set them throughout calculation
// 4) check if node_count_min = 0 for first eigenvalue (either parity) -- done 
// 5) are 'preliminary' and 'precise' types of eigenvalue necessary? clean them? 
// 6) a problem with finding classical turning points when energy is less than minimum of potential -- done

void fixNodeCount( std::vector<Eigenvalue> & eigs );

std::string parityToString( const Parity parity ); 

std::vector<Eigenvalue> calculate_eigenvalues( const Parity parity, Equations & equations, const int NPoints, const double E_min, const double E_max, const int energy_intervals, const double x_lb, const double x_rb, const int i_match_intervals, const double eps, const double eps_tp );

int main( int argc, char * argv[] )
{
    const int NPoints = 5000;
    
    int J, M, channels;
    if ( argc < 4 ) 
    {
        std::cerr << "Program expects (int) J, (int) M, (int) channels" << std::endl;
        exit( 1 );
    }
    else
    {
        J = atoi( argv[1] );
        M = atoi( argv[2] );
        channels = atoi( argv[3] );
        
        std::cout << "(main) accepted arguments. J: " << J << ", M: " << M << "; channels: " << channels << "; NPoints: " << NPoints << std::endl; 
    }

    std::cout << std::fixed << std::setprecision(13);

    //int odd_channels = 0;
    //int even_channels = 0;
    //for ( int L = M; L < channels; ++L )
    //{
        //if ( L % 2 == 0 )
            //even_channels++;
        //else
            //odd_channels++;
    //}
    //std::cout << "Number of even channels: " << even_channels << std::endl;
    //std::cout << "Number of odd channels: " << odd_channels << std::endl;

    // energy interval to search eigenvalues in
    const double E_min = -150.0 / constants::HTOCM;
    const double E_max = -100.0 / constants::HTOCM;
    const double eps = 1.0e-1; // relative precision of preliminary eigenvalue
       
    // interval to search turning points in  
    const double x_lb = 5.0;
    const double x_rb = 32.0;
    const double eps_tp = 1.0e-3; // absolute precision of turning point
    
    int energy_intervals = 10; // number of energy intervals [lg(-E_min), lg(-E_max)] to be divided into
    int i_match_intervals = 20; // number of dots condisered to be i_match

    std::vector<Eigenvalue> tmp, preigs;   
      
    Equations equations( channels, NPoints );
    equations.setAngularMomentum( J, M );
    
    equations.fill_V( 18.0, Parity::EMPTY );
    
    const double BFCT = 16.857630;
    const double URED = 18.928832; // ? is it exact value?
    const double RMLMDA = URED / BFCT;
    
    //std::cout << "R: " << 18.0 << "; U: " << std::endl << equations.get_V() << std::endl;
    
    tmp = calculate_eigenvalues( Parity::EMPTY, equations, NPoints, E_min, E_max, energy_intervals, x_lb, x_rb, i_match_intervals, eps, eps_tp );
    if ( tmp.size() != 0 && tmp[0].get_node_count_min() != 0 )
    {
        std::cerr << "Node count min for first eigenvalue is not zero! Exiting..." << std::endl;
        exit( 1 );
    }

    preigs.insert( preigs.end(), tmp.begin(), tmp.end() );

    //if ( even_channels != 0 )
    //{
        //Equations equations( even_channels, NPoints );
        //equations.setAngularMomentum( J, M );		

        //tmp = calculate_eigenvalues( Parity::EVEN, equations, NPoints, E_min, E_max, energy_intervals, x_lb, x_rb, i_match_intervals, eps, eps_tp );
        //if ( tmp.size() != 0 && tmp[0].get_node_count_min() != 0 )   
        //{ 
            //std::cerr << "Node count min for first eigenvalue is not zero! Exiting..." << std::endl;
            //exit( 1 );
        //} 
        //preigs.insert( preigs.end(), tmp.begin(), tmp.end() ); 
    //}

    //if ( odd_channels != 0 )
    //{
        //Equations equations( odd_channels, NPoints );
        //equations.setAngularMomentum( J, M );

        //tmp = calculate_eigenvalues( Parity::ODD, equations, NPoints, E_min, E_max, energy_intervals, x_lb, x_rb, i_match_intervals, eps, eps_tp );
        //if ( tmp.size() != 0 && tmp[0].get_node_count_min() != 0 )   
        //{ 
            //std::cerr << "Node count min for first eigenvalue is not zero! Exiting..." << std::endl;
            //exit( 1 );
        //} 
        //preigs.insert( preigs.end(), tmp.begin(), tmp.end() ); 
    //}
    
    // сбросить node count для чисел, чтобы отфильтровать их по значению  
    //std::for_each( preigs.begin(), preigs.end(), []( Eigenvalue & e ) { e.dump_node_count(); });
    //std::sort( preigs.begin(), preigs.end() );
    // выставить значения node count по возрастанию
    //fixNodeCount( preigs );

    std::cout << "List of sorted eigenvalues: " << std::endl;
    for ( size_t k = 0; k < preigs.size(); ++k )
        std::cout << preigs[k] << std::endl; 

    const std::string database = "./database/dunker_J" + std::to_string(J) + "M" + std::to_string(M) + "_ch" + std::to_string(channels) + ".db";
    std::cout << "Database filename: " << database << std::endl;

    EigenvalueDatabase EigenDB( database );    
    EigenDB.read();
   
    std::cout << "Database values:" << std::endl;
    EigenDB.show();

    EigenDB.insert( preigs );
    EigenDB.write();    

	return 0;
}


    //EigenfunctionFinder eff( &equations, channels, NPoints );

    //std::vector<Eigen::VectorXd> eigenfunction( NPoints );
    //for ( size_t k = 0; k < NPoints; ++k )
        //eigenfunction[k] = Eigen::VectorXd::Zero( channels );

    //for ( size_t k = 0; k < eigs.size(); ++k )
    //{
        //std::string filename = std::to_string(k) + ".txt";
        //std::cout << "filename: " << filename << std::endl;

        //eff.calculate_eigenfunction( eigenfunction, precise_eigenvalues[k].value, a, b, h, i_match ); 

        //std::ofstream file( filename );
        //file << std::fixed << std::setprecision(14);    

        //double x = a;
        //for ( int k = 0; k < NPoints; ++k, x += h )
            //file << x << " " <<  eigenfunction[k](0) << " " << eigenfunction[k](1) << " " << eigenfunction[k](2) << " " << eigenfunction[k](3) << std::endl;
    //}
    


void fixNodeCount( std::vector<Eigenvalue> & eigs )
{
    for ( size_t k = 0; k < eigs.size(); ++k )
    {
        eigs[k].set_node_count_min( k );
        eigs[k].set_node_count_max( k + 1 );
    }
}

std::string parityToString( const Parity parity ) 
{
    if ( parity == Parity::EVEN )
        return "even";
    else if ( parity == Parity::ODD )
        return "odd";
    else
       return "general"; 
}


std::vector<Eigenvalue> calculate_eigenvalues( const Parity parity, Equations & equations, const int NPoints, const double E_min, const double E_max, const int energy_intervals, const double x_lb, const double x_rb, const int i_match_intervals, const double eps, const double eps_tp )
{
    std::vector<Eigenvalue> preigs;

    std::map<double, std::pair<double, double>> energy_dict = equations.create_energy_dict( parity, E_min, E_max, energy_intervals, x_lb, x_rb, eps_tp );
    bool all_tp_zero = true;
    for ( auto it = energy_dict.begin(); it != energy_dict.end(); ++it )
    {
        std::cout << it->first << " " << it->second.first << " " << it->second.second << std::endl;
        if ( (it->second.first != 0.0) && (it->second.second != 0.0) )
        {
            all_tp_zero = false;
        }
    }
    if ( !all_tp_zero )
    {
        PreliminaryEigenvalueFinder pef( &equations, energy_dict, parity );

        std::vector<Eigenvalue> eigs = pef.findEigenvalues( E_min, E_max );

        std::cout << "Preliminary " + parityToString(parity) + " eigenvalues: " << std::endl;
        for ( size_t k = 0; k < eigs.size(); ++k )
        {
            eigs[k] = pef.convergeToEigenvalue( eigs[k], std::abs(eigs[k].get_value()) * eps );
            std::cout << eigs[k] << std::endl;
        }

        PreciseEigenvalueFinder pref( &equations, energy_dict, parity );
        pref.reserve_space_for_search( i_match_intervals );

        int i_match;
        std::cout << "Precise " + parityToString(parity) + " eigenvalues: " << std::endl;
        for ( size_t k = 0; k < eigs.size(); ++k )
        {
            i_match = pref.effecient_i_match( eigs[k].get_min(), eigs[k].get_max(), NPoints, i_match_intervals );
            std::cout << "(main) i_match: " << i_match << std::endl;
            if ( i_match < 0 )
                continue;

            double eig = pref.precise_eigenvalue_calculation( i_match, eigs[k].get_min(), eigs[k].get_max() );

            preigs.push_back( eigs[k] );
            preigs.back().set_value( eig );

            std::cout << preigs.back() << std::endl; 
        }
    }
    else
    {
        std::cerr << "All turning points are zero! There are no bound levels of " + parityToString(parity) + " parity! Check channels." << std::endl;
    } 

    return preigs;
}


