#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <functional>

#include <Eigen/Dense>
#include <gsl/gsl_sf_coupling.h>

class Functor 
{
public:
    Functor( std::function<int(double)> f ) : f(f)
    {
        calls = 0;
    }

    int operator()( double x )
    {
        ++calls;
        return f( x );
    }

    int get_calls()
    {
        return calls;
    }

    void zero_out()
    {
        calls = 0;
    }
private:
    int calls;
    std::function<int(double)> f;
};

class Solver
{
public:
	Solver( int channels, int NPoints ) : 
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

	void propagateDetR( const double Energy, const double a, const double b, const double h, std::vector<double> & nodesPos )
	{
		// clearing nodes vector!
		nodesPos.clear();

		Rinv = Eigen::MatrixXd::Zero( channels, channels );
		
        double hh12 = h * h / 12.0;

		double detR = std::numeric_limits<double>::max();

		double x = a + h;
		for ( int i = 1; i < NPoints - 1; i++, x += h )
		{
			fill_V( x );
			
			Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
			U = 12.0 * Wmat.inverse() - 10.0 * I;
			R = U - Rinv;
			Rinv = R.inverse();

			detR = R.determinant();
			if ( detR < 0 )
				nodesPos.push_back( x );
		}
	}

    inline int delta( int k, int l ) { return k == l; }

private:
	
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

double bisection( Functor & f, double a, double b, const double eps )
{
    int fa = f( a );
    int fb = f( b );

    if ( fa == fb )
    {
        std::cerr << "Given interval doesn't bracket eigenvalue!" << std::endl;
        exit( 1 );
    }

    double x;
    int fx;

    while( std::abs(a - b) > eps )
    {
        x = 0.5 * (a + b);
        fx = f(x);

        if ( fx == fa )
        {
            a = x;
            fa = fx;
        }
        else
        {
            b = x;
            fb = fx;
        }
    }

    return 0.5 * (a + b);
}

double trisection( Functor & f, double a, double b, const double eps )
{
    int fa = f( a );
    int fb = f( b );

    if ( fa == fb )
    {
        std::cerr << "Given interval doesn't bracket eigenvalue!" << std::endl;
        exit( 1 );
    }

    int fx1, fx2;
    double x1, x2;

    while ( std::abs(b - a) > eps )
    {
        x1 = a + (b - a) / 3.0;
        fx1 = f( x1 );

        if ( fa != fx1 )
        {
            b = x1;
            fb = fx1;
        }
        else
        {
            x2 = b - (b - a) / 3.0;
            fx2 = f( x2 );

            if ( fx1 != fx2 )
            {
                a = x1;
                fa = fx1;
                b = x2;
                fb = fx2;
            }
            else
            {
                a = x2;
                fa = fx2;
            }
        }
    }
    
    return 0.5 * (a + b);
}

double Interpolate2( double x1, double x2, double fx1, double fx2 )
// inverse linear interpolation 
{
    return (x1 * fx2 - x2 * fx1) / (fx2 - fx1);
}

double trisection_plus( Functor & f, double a, double b, const double eps )
{
    int fa = f( a );
    int fb = f( b );

    if ( fa == fb )
    {
        std::cerr << "Given interval doesn't bracket eigenvalue!" << std::endl;
        exit( 1 );
    }

    double lastA, lastB;
    double x1, x2, x3;
    int fx1, fx2, fx3;

    do 
    {
        lastA = a;
        lastB = b;

        x1 = a + (b - a) / 3.0;
        fx1 = ( x1 );

        // case 1: [a, x1] has the root
        if ( fx1 != fa )
        {
            x3 = Interpolate2( a, x1, fa, fx1 );
            fx3 = f( x3 );

            if ( fa != fx3 )
            {
                b = x3;
                fb = fx3;
            }
            else
            {
                a = x3;
                fa = fx3;
                b = x1;
                fb = fx1;
            }
        }
        else
        {
            x2 = a + 2.0 * (b - a) / 3.0;
            fx2 = f( x2 );

            // case 2: [x1, x2] has root
            if ( fx1 != fx2 )
            {
                x3 = Interpolate2( x1, x2, fx1, fx2 );
                fx3 = f( x3 );

                if ( fx1 != fx3 )
                {
                    a = x1;
                    fa = fx1;
                    b = x3;
                    fb = fx3;
                }
                else
                {
                    a = x3;
                    fa = fx3;
                    b = x2;
                    fb = fx2;
                }
            }
            else
            {
                // case 2: [x2, b] has the root
                x3 = Interpolate2( x2, b, fx2, fb );
                fx3 = f( x3 );

                if ( fx2 != fx3 )
                {
                    a = x2;
                    fa = fx2;
                    b = x3;
                    fb = fx3;
                }
                else
                {
                    a = x3;
                    fa = fx3;
                }
            }
        }
        
        if ( (lastA != a) && (std::abs(a - lastA) < eps) )
            break;
        if ( (lastB != b) && (std::abs(b - lastB) < eps) )
            break;
    }
    while ( std::abs(a - b) > eps );

    return a;
}


int main()
{
    const int channels = 4;
    const int J = 7;
    const int M = 0;

    const double a = 6.0;
    const double b = 14.0;
    const int NPoints = 500;
    const double h = (b - a) / (NPoints - 1);

    Solver solver( channels, NPoints );
    solver.setAngularMomentum( J, M );

    std::function<int(double)> f = [&]( double E )
    {
        std::vector<double> nodesPos;
        solver.propagateDetR( E, a, b, h, nodesPos );
        return nodesPos.size(); 
    };
    
    Functor functor( f );

    double E1 = -0.00055;
    double E2 = -0.00054;
    double eps = 1.0e-8;
   
    std::cout << std::fixed << std::setprecision(13);
     
    double E_bis = bisection( functor, E1, E2, eps );
    std::cout << "Bisection energy: " << E_bis << std::endl;
    std::cout << "Function calls: " << functor.get_calls() << std::endl;

    functor.zero_out();

    double E_tri = trisection( functor, E1, E2, eps );
    std::cout << "Trisection energy: " << E_tri << std::endl;
    std::cout << "Function calls: " << functor.get_calls() << std::endl;

    functor.zero_out();

    double E_tri_plus = trisection_plus( functor, E1, E2, eps );
    std::cout << "Trisection-plus energy: " << E_tri_plus << std::endl;
    std::cout << "Function calls: " << functor.get_calls() << std::endl;

    return 0;
}
