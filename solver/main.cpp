#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <map>

#include "./parity.hpp"
#include "./equations.hpp"
#include "./preliminary_eigenvalue_finder.hpp"
#include "./precise_eigenvalue_finder.hpp"
//#include "./eigenfunction_finder.hpp"

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

int main()
{
    std::cout << std::fixed << std::setprecision(13);
	
    const int NPoints = 5000;
    const int J = 7;
    const int M = 0;
    
    //const int channels = J - std::abs(M) + 1;
    int channels = 4;

    Equations equations( channels, NPoints );
    equations.setAngularMomentum( J, M );		

    double E_min = -3.0e-3;
    double E_max = -1.0e-10;
    double eps = 1.0e-1; // relative precision of preliminary eigenvalue
    const double eps_tp = 1.0e-3; // absolute precision of turning point
        
    Parity parity;
    int i_match;
    double eig;

    double x_lb = 5.0;
    double x_rb = 20.0;
    
    bool all_tp_zero;
    int energy_intervals = 10; // number of energy intervals [lg(-E_min), lg(-E_max)] to be divided into
    int i_match_intervals = 100; // number of dots condisered to be i_match
    
    // --------------------------------------------------------------------------------------------------------------
    std::vector<Eigenvalue> even_preigs;
    
    std::map<double, std::pair<double, double>> even_energy_dict = equations.create_energy_dict( Parity::EVEN, E_min, E_max, energy_intervals, x_lb, x_rb, eps_tp );
    std::cout << "Even energy_dict: " << std::endl;
    
    all_tp_zero = true;
    for ( auto it = even_energy_dict.begin(); it != even_energy_dict.end(); ++it )
    {
        std::cout << it->first << " " << it->second.first << " " << it->second.second << std::endl;
        if ( (it->second.first != 0.0) && (it->second.second != 0.0) )
        {
           all_tp_zero = false;
        }
    }
    if ( !all_tp_zero )
    {
        parity = Parity::EVEN; 
        PreliminaryEigenvalueFinder even_pef( &equations, even_energy_dict, parity );

        std::vector<Eigenvalue> even_eigs = even_pef.findEigenvalues( E_min, E_max );
        std::cout << "Preliminary even eigenvalues: " << std::endl;
        for ( size_t k = 0; k < even_eigs.size(); ++k )
        {
            even_eigs[k] = even_pef.convergeToEigenvalue( even_eigs[k], std::abs(even_eigs[k].get_value()) * eps );
            std::cout << even_eigs[k] << std::endl;
        }

        PreciseEigenvalueFinder even_pref( &equations, even_energy_dict, parity );
        even_pref.reserve_space_for_search( i_match_intervals );

        std::cout << "Precise even eigenvalues: " << std::endl;
        for ( size_t k = 0; k < even_eigs.size(); ++k )
        {
            i_match = even_pref.effecient_i_match( even_eigs[k].get_min(), even_eigs[k].get_max(), NPoints, i_match_intervals );
            std::cout << "(main) i_match: " << i_match << std::endl;
            if ( i_match < 0 )
                continue;

            eig = even_pref.precise_eigenvalue_calculation( i_match, even_eigs[k].get_min(), even_eigs[k].get_max() );

            even_preigs.push_back( even_eigs[k] );
            even_preigs.back().set_value( eig );

            std::cout << even_preigs.back() << std::endl; 
        }
    }
    else
    { 
        std::cerr << "All turning points are zero! There are no bound levels of EVEN parity! Check channels." << std::endl;
    }
    // --------------------------------------------------------------------------------------------------------------

    // --------------------------------------------------------------------------------------------------------------
    std::vector<Eigenvalue> odd_preigs;
    
    std::map<double, std::pair<double, double>> odd_energy_dict = equations.create_energy_dict( Parity::ODD, E_min, E_max, energy_intervals, x_lb, x_rb, eps_tp );
    all_tp_zero = true;
    std::cout << "Odd energy_dict: " << std::endl;
    for ( auto it = odd_energy_dict.begin(); it != odd_energy_dict.end(); ++it )
    {
        std::cout << it->first << " " << it->second.first << " " << it->second.second << std::endl; 
        if ( (it->second.first != 0.0) && (it->second.second != 0.0) )
        {
           all_tp_zero = false;
        }
    }
    
    if ( !all_tp_zero )
    {
        parity = Parity::ODD;
        PreliminaryEigenvalueFinder odd_pef( &equations, odd_energy_dict, parity );

        std::vector<Eigenvalue> odd_eigs = odd_pef.findEigenvalues( E_min, E_max );
        std::cout << "Preliminary odd eiganvalues: " << std::endl;
        for ( size_t k = 0; k < odd_eigs.size(); ++k )
        {
            odd_eigs[k] = odd_pef.convergeToEigenvalue( odd_eigs[k], std::abs(odd_eigs[k].get_value()) * eps );
            std::cout << odd_eigs[k] << std::endl;
        }

        PreciseEigenvalueFinder odd_pref( &equations, odd_energy_dict, parity );
        odd_pref.reserve_space_for_search( i_match_intervals );

        std::cout << "Precise odd eigenvalues: " << std::endl;
        for ( size_t k = 0; k < odd_eigs.size(); ++k )
        {
            i_match = odd_pref.effecient_i_match( odd_eigs[k].get_min(), odd_eigs[k].get_max(), NPoints, i_match_intervals );
            std::cout << "(main) i_match: " << i_match << std::endl;
            if ( i_match < 0 )
                continue;

            eig = odd_pref.precise_eigenvalue_calculation( i_match, odd_eigs[k].get_min(), odd_eigs[k].get_max() );

            odd_preigs.push_back( odd_eigs[k] );
            odd_preigs.back().set_value( eig );

            std::cout << odd_preigs.back() << std::endl;
        }
    }
    else
    {
        std::cerr << "All turning points are zero! There are no bound levels of ODD parity! Check channels." << std::endl;
    } 
    // --------------------------------------------------------------------------------------------------------------    

    if ( (even_preigs.size() != 0 && even_preigs[0].get_node_count_min() != 0) || (odd_preigs.size() != 0.0 && odd_preigs[0].get_node_count_min() != 0) )
    {
        std::cout << "Node count min for first eigenvalue is not zero!" << std::endl << "Exiting without saving to database..." << std::endl;
        exit( 1 );
    }
    
    std::vector<Eigenvalue> preigs; // full list of eigenvalues
    preigs.insert( preigs.end(), even_preigs.begin(), even_preigs.end() );
    preigs.insert( preigs.end(), odd_preigs.begin(), odd_preigs.end() );
    
    // сбросить node count для чисел, чтобы отфильтровать их по значению  
    std::for_each( preigs.begin(), preigs.end(), []( Eigenvalue & e ) { e.dump_node_count(); });
    std::sort( preigs.begin(), preigs.end() );
    // выставить значения node count по возрастанию
    fixNodeCount( preigs );

    std::cout << "List of sorted eigenvalues: " << std::endl;
    for ( size_t k = 0; k < preigs.size(); ++k )
        std::cout << preigs[k] << std::endl; 

    const std::string database = "./database/J" + std::to_string(J) + "M" + std::to_string(M) + ".db";
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

