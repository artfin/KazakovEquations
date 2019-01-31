#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <map>
#include "./Ar-HCl_pes_coefs.hpp"
//#include "./Hutson_Ar_HCl_pes_coefs.hpp"

#include <Eigen/Dense>

#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include "./constants.hpp"
#include "./parity.hpp"

class Equations
{
public:
	explicit Equations( int channels, int NPoints ); 

	void setAngularMomentum( const int J, const int M );
    void init_matrices();
    void reset_channels( const int channels );

	void fill_W_elements( const double R );
    std::vector<double> get_W() const { return W; }     

	double angularMatrixElements( const int P, const int l, const int Lprime );
	
	double compute_W_sum( const int P, const int Lprime);

	void fill_V( const double r, const Parity parity );

	Eigen::MatrixXd get_V( ) { return V; }
	int delta( int k, int l ) { return k == l; }
    int get_channels( ) const { return channels; }
    int get_M() const { return M; }
    int get_J() const { return J; }

    void propagateDetR( const Parity parity, const double Energy, const double a, const double b, const double h, std::vector<double> & nodesPos );
    void propagateForward( const Parity parity, const double Energy, const double a, const double h, const int i_match, Eigen::MatrixXd & resRm, bool save = false );
    void propagateForwardFull( const Parity parity, const double Energy, const double a, const double h, std::vector<Eigen::MatrixXd> & Rm_vector, const int step );
    void propagateBackward( const Parity parity, const double Energy, const double b, const double h, const int i_match, Eigen::MatrixXd & resRmp1, bool save = false ); 
    void propagateBackwardFull( const Parity parity, const double Energy, const double b, const double h, std::vector<Eigen::MatrixXd> & Rmp1_vector, const int step ); 

    int countEigenvalues( const Parity parity, const double Energy, const double a, const double b, const double h );  
    int countNegativeDiagonalElements( );

    double brent( std::function<double(double)> f, double xb1, double xb2, const double eps );
    Eigen::VectorXd adiabatic_potential_vector( const Parity parity, const double r );
    double adiabatic_potential_component( const Parity parity, const double r, const int k );
    std::pair<double, double> find_classical_turning_points( const Parity parity, const double Energy, double x_lb, double x_rb, const double eps );
    std::map<double, std::pair<double, double>> create_energy_dict( const Parity parity, const double E_min, const double E_max, const int energy_intervals, const double x_lb, const double x_rb, const double eps );
    std::pair<double, double> interpolate( const double e, const std::map<double, std::pair<double, double>> & energy_dict );
    void calculate_boundaries( std::pair<double, double> const & tp, double * a, double * b, double * h );

    std::vector<Eigen::MatrixXd> Rinv_vector;
    std::vector<Eigen::MatrixXd> Winv_vector;

    Eigen::MatrixXd const& get_V() const { return V; }

private:
	const double hbar = 1.0;
    //const double mu = 38183.0;
    //const double Inten = 14579.0 * std::pow(4.398, 2.0);
    
    const double RAMTOAMU = 1822.888485332;
    const double H1_RAM = 1.00782503223;
    const double AR40_RAM = 39.9623831237;
    const double CL35_RAM = 34.968852682;

    const double HCL_MU = H1_RAM * CL35_RAM / (H1_RAM + CL35_RAM) * RAMTOAMU;

    //const double mu = AR40_RAM * (H1_RAM + CL35_RAM) / (AR40_RAM + H1_RAM + CL35_RAM) * RAMTOAMU;
    const double mu = 34505.15; 

    // atomic length unit to meters == a0
    const double ALU = 0.52917721067 * 1E-10;

    const double HCL_LEN = 0.127e-9 / ALU; 
    //const double Inten = HCL_MU * std::pow(HCL_LEN, 2.0);
    //
    const double alpha_inv = 137.035999139;
    const double B = 10.44019 * 0.529177e-8; // a.l.u
    const double Inten = 1.0 / (4.0 * M_PI * B * alpha_inv); 

	int channels;
    int NPoints;

	int M = 0;
    int M_abs = 0;
    int J = -1;

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
