#include <iostream>
#include <iomanip>
#include <fstream>

#include "./equations.hpp"
#include "./preliminary_eigenvalue_finder.hpp"
#include "./precise_eigenvalue_finder.hpp"
#include "./eigenfunction_finder.hpp"

int main()
{
	const int NPoints = 500;
	const int channels = 4;
	Equations equations( channels, NPoints );

	const int J = 7;
	const int M = 0;
	equations.setAngularMomentum( J, M );		

	PreliminaryEigenvalueFinder pef( &equations, NPoints );
   
    double a = 5.0;
    double b = 50.0;
    double h = (b - a) / (NPoints - 1);

    //int nodes = equations.countEigenvalues( -1.0e-12, a, b, h ); 
    //std::cout << "Eigenvalues: " << nodes << std::endl; 

    //double E_min = -5.5e-4;
    //double E_max = -1.0e-10;
    //double eps = 1.0e-8;

    //std::cout << std::fixed << std::setprecision(13);

    //std::vector<Eigenvalue> eigs = pef.findEigenvalues( E_min, E_max, a, b, h );
    //std::cout << "-----------------------------------------" << std::endl;
    //std::cout << "Found eigenvalues: " << std::endl;
    //for ( size_t k = 0; k < eigs.size(); ++k )
        //std::cout << eigs[k] << std::endl;

    
    //std::cout << "-----------------------------------------" << std::endl;
    //std::cout << "Converged eigenvalues: " << std::endl;
    //for ( size_t k = 0; k < eigs.size(); ++k )
    //{
        //eigs[k] = pef.convergeToEigenvalue( eigs[k], a, b, h, eps );
        //std::cout << eigs[k] << std::endl;
    //}

    //PreciseEigenvalueFinder preciseEF( &equations, NPoints );

    //int i_match;
    //double eig;

    //std::vector<Eigenvalue> precise_eigenvalues;
    //for ( size_t k = 0; k < eigs.size(); ++k )
    //{
        //i_match = preciseEF.find_i_match( a, b, h, eigs[k].min, eigs[k].max );
        //eig = preciseEF.precise_eigenvalue_calculation( a, b, h, i_match, eigs[k].min, eigs[k].max );

        //precise_eigenvalues.emplace_back( eig, eigs[k].node_count_min, eigs[k].node_count_max );
    //}

    //std::cout << "----------------------------------------------" << std::endl;
    //std::cout << "List of precise eigenvalues: " << std::endl;
    //for ( size_t k = 0; k < precise_eigenvalues.size(); ++k )
        //std::cout << precise_eigenvalues[k] << std::endl;
    //std::cout << "----------------------------------------------" << std::endl;

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


	return 0;
}
