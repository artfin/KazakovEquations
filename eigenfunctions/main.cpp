#include <iostream>
#include <vector>
#include <iomanip>

#include <Eigen/Dense>

#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

class EigenfunctionSolver
{
public:
    EigenfunctionSolver( int channels, int NPoints ) :
        channels(channels), NPoints(NPoints)
    {
        W.resize(6);

		I = Eigen::MatrixXd::Identity(channels, channels);
		
		V = Eigen::MatrixXd::Zero( channels, channels ); 
		R = Eigen::MatrixXd::Zero(channels, channels);
		Rinv = Eigen::MatrixXd::Zero(channels, channels);
		Wmat = Eigen::MatrixXd::Zero(channels, channels);
		U = Eigen::MatrixXd::Zero(channels, channels);
    }

	void setAngularMomentum( const int J, const int M )
	{
		this->J = J;
		this->M = M;
	}

	void fill_W_elements( const double R )
	{
		W[0] = 22.4247 * std::exp(-0.716288*R-0.0869136*R*R);
		
		if (R >= 6.32925)
			W[0] -= 114.5 / std::pow(R, 6) + 2380.0 / std::pow(R, 8);
		else
			W[0] += (-0.599056 * std::exp(-0.650479*R-0.0320299*R*R));
			
		W[1] = 63.5744 * std::exp(-0.811806*R-0.075313*R*R);

		if (R < 6.90026)
			W[1] -= 0.207391 * std::exp(-0.620877*R-0.0310717*R*R);
		else
			W[1] -= 26.6 / std::pow(R, 6) + 2080.0 / std::pow(R,8);

		W[2] = 99.1128 * std::exp(-1.17577*R-0.0477138*R*R);

		if (R < 6.96450)
			W[2] -= 0.0523497 * std::exp(-0.73534*R-0.0296750*R*R);
		else 
			W[2] -= 410 / std::pow(R,8);

		W[3] = 318.652 * std::exp(-1.88135*R) - 0.0374994 * std::exp(-1.05547*R-0.0182219*R*R);

		W[4] = 332.826 * std::exp(-2.14596*R) - 0.0137376 * std::exp(-1.03942*R-0.0494781*R*R);

		W[5] = 435.837 * std::exp(-2.44616*R) - 0.108283 * std::exp(-2.04765*R);
	}
	
	double angularMatrixElements( const int P, const int l, const int Lprime )
	// gsl_sf_coupling_3j принимает удвоенные значения спинов и проекций спинов
	// возвращает значение матричного элемента <PM|P_l|L'M>
	{
		return std::pow(-1.0, M) * std::sqrt( (2.0*P+1.0) * (2.0*Lprime+1.0) ) * gsl_sf_coupling_3j(2*Lprime, 2*l, 2*P, 0, 0, 0 ) * gsl_sf_coupling_3j(2*Lprime, 2*l, 2*P, 2*M, 0, -2*M );
	}
	
	double compute_W_sum( const int P, const int Lprime)
	// возвращает значение суммы \sum_l W_l * <PM|P_l|L'M>
	{
		double result = 0.0;
		for ( size_t i = 0; i < W.size(); ++i )
			result += W[i] * angularMatrixElements( P, 2 * i, Lprime );
		
		return result;
	}

	void fill_V( const double r )
	{
		// заполняем элементы вектора W
		fill_W_elements( r );
	
		for ( int P = 0; P < V.rows(); ++P )
			for ( int Lprime = 0; Lprime < V.cols(); ++Lprime )
                //V(P, Lprime) = compute_W_sum( 2*P, 2*Lprime ) + delta( 2*P, 2*Lprime ) * (hbar*hbar*2.0*P*(2.0*P+1.0)/2.0/Inten + hbar*hbar/2.0/mu/r/r*(J*(J+1) + 2.0*P*(2.0*P+1.0)));
                V(P, Lprime) = compute_W_sum( P, Lprime ) + delta(P, Lprime) * (hbar*hbar*P*(P+1)/2.0/Inten + hbar*hbar/2.0/mu/r/r*(J*(J+1) + P*(P+1)));
	}

	Eigen::MatrixXd get_V( ) { return V; }
	
    void propagateForward( const double Energy, const double a, const double h, const int i_match, Eigen::MatrixXd & resRm ) 
	{
		double hh12 = h * h / 12.0;
		double x = a + h;

		Rinv = Eigen::MatrixXd::Zero( channels, channels );

		for ( int i = 1; i < i_match; i++, x += h )
		{
			fill_V( x );

			Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
			U = 12.0 * Wmat.inverse() - 10.0 * I;
			R = U - Rinv;
			Rinv = R.inverse();
		}	
	
		resRm = R;
	}
	
	void propagateBackward( const double Energy, const double b, const double h, const int i_match, Eigen::MatrixXd & resRmp1 ) 
	// resRmp1 = R_{m + 1}
	{
		double hh12 = h * h / 12.0;
		double x = b - h;
		
		Rinv = Eigen::MatrixXd::Zero( channels, channels );

		for ( int i = NPoints - 1; i > i_match; i--, x -= h )
		{
			fill_V( x );

			Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
			U = 12.0 * Wmat.inverse() - 10.0 * I;
			R = U - Rinv;
			Rinv = R.inverse();
		}

		resRmp1 = R;	
	}
	
    double brent( std::function<double(double)> f, double xb1, double xb2 )
	{
		int status;
		int iter = 0, max_iter = 100;
	
		const gsl_root_fsolver_type * T;
		T = gsl_root_fsolver_brent;
	
		gsl_root_fsolver * s;
		s = gsl_root_fsolver_alloc( T );

		gsl_function F = 
		{
			[](double d, void * vf) -> double {
				auto& f = *static_cast<std::function<double(double)>*>(vf);
				return f(d);
			},
			&f
		};

		gsl_root_fsolver_set( s, &F, xb1, xb2 );

		double x_lo, x_hi, r = 0;

		const double eps_abs = 1.0e-13;
		const double eps_rel = 1.0e-13;

		do
		{
			iter++;

			status = gsl_root_fsolver_iterate( s );
			r = gsl_root_fsolver_root( s );

			x_lo = gsl_root_fsolver_x_lower( s );
			x_hi = gsl_root_fsolver_x_upper( s );
			status = gsl_root_test_interval( x_lo, x_hi, eps_abs, eps_rel );

			//if ( status == GSL_SUCCESS )
				//printf( "Converged:\n");

			//printf( "%5d [%.12f, %.12f] %.12f %.12f\n",
				//iter, x_lo, x_hi, r, x_hi - x_lo );	
		}
		while ( status == GSL_CONTINUE && iter < max_iter );	

		gsl_root_fsolver_free( s );

		return r;
	}

	double D( const double a, const double b, const double h, const int i_match, const double E )
	{
		Eigen::MatrixXd Rm, Rmp1;

		propagateForward( E, a, h, i_match, Rm );
		propagateBackward( E, b, h, i_match, Rmp1 );
		
		return (Rm - Rmp1.inverse()).determinant();
	}

    void calculate_eigenfunction( const double a, const double b, const double h, const int i_match, const double E1, const double E2 )
    {
        std::function<double(double)> __D__ = [=]( const double E ) { return D(a, b, h, i_match, E); };
        double eig = brent( __D__, E1, E2 ); 
        std::cout << std::fixed << std::setprecision(14);
        std::cout << "brent eigenvalue: " << eig << std::endl;   
        std::cout << "Det: " << __D__(eig) << std::endl;

        Eigen::MatrixXd Rm, Rmp1;
        propagateForward( eig, a, h, i_match, Rm );
        propagateBackward( eig, b, h, i_match, Rmp1 );
        Eigen::MatrixXd diff = Rm - Rmp1.inverse();
        std::cout << "diff: " << std::endl << diff << std::endl;

        Eigen::FullPivLU<Eigen::MatrixXd> lu( diff );
        lu.setThreshold( 1.0e-8 );

        Eigen::MatrixXd ker = lu.kernel();
        std::cout << "kernel: " << std::endl << ker << std::endl;

        //std::cout << "E1 = " << E1 << "; matrix: " << std::endl << diff1 << std::endl;

        //propagateForward( E2, a, h, i_match, Rm );
        //propagateBackward( E2, b, h, i_match, Rmp1 );

        //diff2 = Rm - Rmp1.inverse();
        //std::cout << "E2 = " << E2 << "; matrix: " << std::endl << diff2 << std::endl;
    
        //std::cout << "difference between matrices: " << std::endl << diff1 - diff2 << std::endl;
    }

private:
	inline int delta( int k, int l ) { return k == l; }
	
	const double hbar = 1.0;
	const double mu = 38183.0;
	const double Inten = 14579.0 * std::pow(4.398, 2.0);

	int channels;
	int NPoints;

	int M;
	int J;

	std::vector<double> W;

	Eigen::MatrixXd V;
	Eigen::VectorXd Veig;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

	Eigen::MatrixXd Wmat, U, R, Rinv;
	Eigen::MatrixXd	I;
};

int main()
{
    int channels = 4;
    int NPoints = 5000;

    const double a = 6.0;
    const double b = 14.0;
    const double h = (b - a) / (NPoints - 1);

    EigenfunctionSolver solver( channels, NPoints ); 
    int J = 7;
    int M = 0;
    solver.setAngularMomentum( J, M );

    //double E1 = -0.000541;
    //double E2 = -0.000540;
    //double E1 = -0.0001618283108;
    //double E2 = -0.0001618283100;
    double E1 = -0.00009708673;
    double E2 = -0.00009708672;
    
    //std::cout << std::fixed << std::setprecision(13);
    //std::cout << "(D) E1 = " << E1 << "; D = " << solver.D(a, b, h, i_match, E1) << std::endl;
    //std::cout << "(D) E2 = " << E2 << "; D = " << solver.D(a, b, h, i_match, E2) << std::endl;

    double D1, D2;
    int i_match;
    const int parts = 10;
    for ( int k = 0; k < parts; ++k )
    {
        i_match = NPoints / parts * k;
        std::cout << "i_match: " << i_match << std::endl;
        
        D1 = solver.D(a, b, h, i_match, E1);
        D2 = solver.D(a, b, h, i_match, E2);

        if ( D1 * D2 < 0.0 )
        {
            std::cout << "D1 * D2 < 0.0! i_match = " << i_match << " is appropriate!" << std::endl;
            break;
        }
    }

    solver.calculate_eigenfunction( a, b, h, i_match, E1, E2 );

    return 0;
}

