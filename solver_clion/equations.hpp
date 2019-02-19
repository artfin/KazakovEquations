#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <map>

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
	explicit Equations( unsigned int channels, unsigned int NPoints );

	void setAngularMomentum( int J, int M );
    void init_matrices();
    void reset_channels( unsigned int channels );

	void fill_W_elements( double R );
    std::vector<double> get_W() const { return W; }     

	double angularMatrixElements( int P, int l, int Lprime );
	
	double compute_W_sum( int P, int Lprime);

	void fill_V( double r, Parity parity );

	Eigen::MatrixXd get_V( ) { return V; }
	int delta( int k, int l ) { return k == l; }
    int get_channels( ) const { return channels; }
    int get_M() const { return M; }
    int get_J() const { return J; }

    void propagateDetR( Parity parity, double Energy, double a, double b, double h, std::vector<double> & nodesPos );
    void propagateForward( Parity parity, double Energy, double a, double h, int i_match, Eigen::MatrixXd & resRm, bool save = false );
    void propagateForwardFull( Parity parity, double Energy, double a, double h, std::vector<Eigen::MatrixXd> & Rm_vector, int step );
    void propagateBackward( Parity parity, double Energy, double b, double h, int i_match, Eigen::MatrixXd & resRmp1, bool save = false );
    void propagateBackwardFull( Parity parity, double Energy, double b, double h, std::vector<Eigen::MatrixXd> & Rmp1_vector, int step );

    int countEigenvalues( Parity parity, double Energy, double a, double b, double h );
    int countNegativeDiagonalElements( );

    double brent( std::function<double(double)> f, double xb1, double xb2, double eps );
    Eigen::VectorXd adiabatic_potential_vector( Parity parity, double r );
    double adiabatic_potential_component( Parity parity, double r, int k );
    std::pair<double, double> find_classical_turning_points( Parity parity, double Energy, double x_lb, double x_rb, double eps );
    std::map<double, std::pair<double, double>> create_energy_dict( Parity parity, double E_min, double E_max, int energy_intervals,
    																double x_lb, double x_rb, double eps );
    std::pair<double, double> interpolate( double e, const std::map<double, std::pair<double, double>> & energy_dict );
    void calculate_boundaries( std::pair<double, double> const & tp, double * a, double * b, double * h );

    std::vector<Eigen::MatrixXd> Rinv_vector;
    std::vector<Eigen::MatrixXd> Winv_vector;

    Eigen::MatrixXd const& get_V() const { return V; }

private:
	// Reduced Planck constant is left for consistency
	const double hbar = 1.0;

   	// CO2-Ar from old files; REDO WITH RAM
    //const double mu = 38183.0;
    //const double Inten = 14579.0 * std::pow(4.398, 2.0);
    
    const double RAMTOAMU = 1822.888485332;

    const double H1_RAM = 1.00782503223;
    const double AR40_RAM = 39.9623831237;
    const double CL35_RAM = 34.968852682;
    const double HCL_MU = H1_RAM * CL35_RAM / (H1_RAM + CL35_RAM) * RAMTOAMU;
    //const double mu = AR40_RAM * (H1_RAM + CL35_RAM) / (AR40_RAM + H1_RAM + CL35_RAM) * RAMTOAMU;

    // taken from Manolopoulos Thesis
    const double mu = 34505.15; 

    const double HCL_LEN = 0.127e-9 / constants::ALU;
    //const double Inten = HCL_MU * std::pow(HCL_LEN, 2.0);
    //
    const double B = 10.44019 * 0.529177e-8; // a.l.u
    const double Inten = 1.0 / (4.0 * M_PI * B * constants::alpha_inv);

	unsigned int channels;
    unsigned int NPoints;

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
