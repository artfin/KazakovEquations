#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <gsl/gsl_sf_coupling.h>

class Equations
{
public:
	explicit Equations( int channels, int NPoints ); 

	void setAngularMomentum( const int J, const int M );
	
	void fill_W_elements( const double R );
	
	double angularMatrixElements( const int P, const int l, const int Lprime );
	
	double compute_W_sum( const int P, const int Lprime);

	void fill_V( const double r );

	Eigen::MatrixXd get_V( ) { return V; }
	int delta( int k, int l ) { return k == l; }

    void propagateDetR( const double Energy, const double a, const double b, const double h, std::vector<double> & nodesPos );
    void propagateForward( const double Energy, const double a, const double h, const int i_match, Eigen::MatrixXd & resRm, bool save = false );
    void propagateBackward( const double Energy, const double b, const double h, const int i_match, Eigen::MatrixXd & resRmp1, bool save = false ); 

    int countEigenvalues( const double Energy, const double a, const double b, const double h );  
    int countNegativeDiagonalElements( );

    std::vector<Eigen::MatrixXd> Rinv_vector;
    std::vector<Eigen::MatrixXd> Winv_vector;

private:
	const double hbar = 1.0;
	const double mu = 38183.0;
	const double Inten = 14579.0 * std::pow(4.398, 2.0);
	
	int channels;
    int NPoints;

	int M;
	int J;

	std::vector<double> W;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    Eigen::VectorXd eigenvalues;

	Eigen::MatrixXd V;
	Eigen::MatrixXd I;
	Eigen::MatrixXd R;
	Eigen::MatrixXd Rinv;
	Eigen::MatrixXd Wmat;
	Eigen::MatrixXd U;

};
