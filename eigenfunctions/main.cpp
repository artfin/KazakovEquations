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

		Rinv_vector.resize( NPoints );
		for ( int k = 0; k < NPoints; ++k )
			Rinv_vector[k] = Eigen::MatrixXd::Zero( channels, channels );

		Winv_vector.resize( NPoints );
		for ( int k = 0; k < NPoints; ++k )
			Winv_vector[k] = Eigen::MatrixXd::Zero( channels, channels );
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
	
    void propagateForward( const double Energy, const double a, const double h, const int i_match, Eigen::MatrixXd & resRm, bool save = false ) 
	{
		double hh12 = h * h / 12.0;
		double x = a + h;

		Rinv = Eigen::MatrixXd::Zero( channels, channels );

		for ( int i = 1; i <= i_match; i++, x += h )
		{
			//std::cout << "(propagateForward) i: " << i << "; x: " << x << std::endl; 
			fill_V( x );

			Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
			U = 12.0 * Wmat.inverse() - 10.0 * I;
			R = U - Rinv;
			Rinv = R.inverse();

			//std::cout << "(propagateForward) Rinv: " << std::endl << Rinv << std::endl;

			if ( save )
			{
				Rinv_vector[i] = Rinv;
				Winv_vector[i] = Wmat.inverse();
			}	
		}
	
		resRm = R;
	}

	void propagateBackward( const double Energy, const double b, const double h, const int i_match, Eigen::MatrixXd & resRmp1, bool save = false ) 
	// resRmp1 = R_{m + 1}
	{
		double hh12 = h * h / 12.0;
		double x = b - h;
		
		Rinv = Eigen::MatrixXd::Zero( channels, channels );

		for ( int i = NPoints - 2; i > i_match; i--, x -= h )
		{
			//std::cout << "(propagateBackward) i: " << i << "; x: " << x << std::endl;
			fill_V( x );
			
			Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
			U = 12.0 * Wmat.inverse() - 10.0 * I;
			R = U - Rinv;
			Rinv = R.inverse();

			//std::cout << "(propagateBackward) Rinv: " << std::endl << Rinv << std::endl;

			if ( save )
			{
				Rinv_vector[i] = Rinv;
				Winv_vector[i] = Wmat.inverse();
			}
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
	
        std::cout << std::setprecision(14) << std::fixed;
        std::cout << "(D) E: " << E << std::endl;
        std::cout << "(D) Rm: " << std::endl << Rm << std::endl;
            
		double v = (Rm - Rmp1.inverse()).determinant();
	    std::cout << "(D) value: " << v << std::endl;

        return v;
    }

	double precise_eigenvalue_calculation( const double a, const double b, const double h, const int i_match, const double E1, const double E2 )
	{
        std::function<double(double)> __D__ = [=]( const double E ) { return D(a, b, h, i_match, E); };
        return brent( __D__, E1, E2 ); 
	}

	std::vector<Eigen::VectorXd> calculate_eigenfunction( const double E, const double a, const double b, const double h, int i_match  )
    {
		//i_match = 3;
		//std::cout << "i_match: " << i_match << std::endl;

        std::cout << std::fixed << std::setprecision(14);
        std::cerr << std::fixed << std::setprecision(14);
		std::cerr << "(calculate_eigenfunction) Eigenvalue: " << E << std::endl;   

        Eigen::MatrixXd Rm, Rmp1;
        propagateForward( E, a, h, i_match, Rm, true );
        propagateBackward( E, b, h, i_match, Rmp1, true );
		Eigen::MatrixXd diff = Rm - Rmp1.inverse();
		
		//std::cout << "diff: " << std::endl << diff << std::endl << std::endl;

        Eigen::FullPivLU<Eigen::MatrixXd> lu( diff );
        lu = lu.setThreshold( 1.0e-5 );

        Eigen::MatrixXd ker = lu.kernel();
		//std::cout << "kernel: " << std::endl << ker << std::endl;
		
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
			//std::cout << "k: " << k << "; Rinv: " << std::endl << Rinv_vector[k] << std::endl;
			fn = Rinv_vector[k] * fn_next;
			psin = Winv_vector[k] * fn;
			eigfunc[k] = psin;	

			//std::cout << "eigfunc[" << k << "]: " << std::endl << eigfunc[k] << std::endl;	
		}	

		for ( int k = i_match + 1; k < NPoints - 1; ++k, fn_prev = fn )
		{
			//std::cout << "k: " << k << "; Rinv: " << std::endl << Rinv_vector[k] << std::endl;
			fn = Rinv_vector[k] * fn_prev;
			psin = Winv_vector[k] * fn;
			eigfunc[k] = psin;	
			
			//std::cout << "eigfunc[" << k << "]: " << std::endl << eigfunc[k] << std::endl;	
		}
		
		psin = Winv_vector[i_match] * fn_matching_point;
		eigfunc[i_match] = psin;

		//std::cout << "eigfunc[" << i_match << "]: " << std::endl << eigfunc[i_match] << std::endl;

		return eigfunc;	
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

	std::vector<Eigen::MatrixXd> Rinv_vector;
	std::vector<Eigen::MatrixXd> Winv_vector;
};

int main()
{
    int channels = 4;
    int NPoints = 500;

    const double a = 6.0;
    const double b = 14.0;
    const double h = (b - a) / (NPoints - 1);

    EigenfunctionSolver solver( channels, NPoints ); 
    int J = 7;
    int M = 0;
    solver.setAngularMomentum( J, M );

	// first eigenvalue: 0 nodes
    double E1 = -0.000541;
    double E2 = -0.000540;
	// second eigenvalue: 1 node
	//double E1 = -0.000405;
	//double E2 = -0.000404;
	// third eigenvalue: 0 nodes
	//double E1 = -0.00032;
	//double E2 = -0.00031;
	// fourth eigenvalue: 
	//double E1 = -0.000291;
	//double E2 = -0.000290;
   	// fifth eigenvalue:
	//double E1 = -0.0002182;
	//double E2 = -0.0002181;
	// 6th eigenvalue
	//double E1 = -0.000200;
	//double E2 = -0.000199;

    double D1, D2;
    int i_match;
    const int parts = 10;
    for ( int k = 0; k < parts; ++k )
    {
        i_match = NPoints / parts * k;
		//std::cout << "i_match: " << i_match << std::endl;
        
        D1 = solver.D(a, b, h, i_match, E1);
        D2 = solver.D(a, b, h, i_match, E2);

        if ( D1 * D2 < 0.0 )
        {
			std::cerr << "D1 * D2 < 0.0! i_match = " << i_match << " is appropriate!" << std::endl;
            break;
        }
    }

	double eig = solver.precise_eigenvalue_calculation( a, b, h, i_match, E1, E2 );
	//double eig = -0.00009708672407; 
	std::vector<Eigen::VectorXd> eigfunc = solver.calculate_eigenfunction( eig, a, b, h, i_match );

	double x = a;
	for ( int k = 0; k < NPoints; ++k, x += h )
		std::cout << x << " " <<  eigfunc[k](0) << " " << eigfunc[k](1) << " " << eigfunc[k](2) << " " << eigfunc[k](3) << std::endl;

    return 0;
}

