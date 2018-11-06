#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <functional>

#include <Eigen/Dense>
//#include <Eigen/src/Core/util/DisableStupidWarnings.h>

#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

struct Eigenvalue
{
	explicit Eigenvalue()
	{
	}

    Eigenvalue( double min, double max, int node_count_min, int node_count_max, double value ) :
        min(min), max(max), node_count_min(node_count_min), node_count_max(node_count_max), value(value)
    {
        degeneracy = node_count_max - node_count_min;
    }

	Eigenvalue( double min, double max, int node_count_min, int node_count_max ) :
		min(min), max(max), node_count_min(node_count_min), node_count_max(node_count_max) 
	{
		value = 0.5 * (min + max);
		degeneracy = node_count_max - node_count_min;
	}

    double min;
    double max;

    int node_count_min;
    int node_count_max;

    double value;
    int degeneracy;
};

std::ostream& operator<< (std::ostream &o, const Eigenvalue & eig)
{
	o << "(Eigenvalue) min: " << eig.min << "; max: " << eig.max << "; value: " << eig.value << "; degeneracy: " << eig.degeneracy << std::endl;
  	return o;
}

struct Interval
{
    Interval( double min, double max, int node_count_min, int node_count_max ) : min(min), max(max), node_count_min(node_count_min), node_count_max(node_count_max) { }

    double min; // нижняя граница энергетического интервала
    double max; // верхняя граница энергетического интервала

    int node_count_min; // количество узлов функции |R|(x) при E = min
    int node_count_max; // количество узлов функции |R|(x) при E = max
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

	Eigen::MatrixXd get_V( ) { return V; }
	
	void propagateDetR_testing( const double Energy, std::vector<double> & xvec, std::vector<double> & detRvec, const double a, const double b, const double h ) 
	{
		Rinv = Eigen::MatrixXd::Zero( channels, channels );
		
        double hh12 = h * h / 12.0;
			
		xvec[0] = a;
		xvec[NPoints - 1] = b;
		detRvec[0] = std::numeric_limits<double>::max();
		detRvec[NPoints - 1] = std::numeric_limits<double>::max();

		double x = a + h;
		for ( int i = 1; i < NPoints - 1; i++, x += h )
		{
			fill_V( x );

			Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
			U = 12.0 * Wmat.inverse() - 10.0 * I;
			R = U - Rinv;
			Rinv = R.inverse();

			xvec[i] =  x;
			detRvec[i] = R.determinant();	
		}	
	}

	void propagateDetR_backward_testing( const double Energy, std::vector<double> & xvec, std::vector<double> & detRvec, const double a, const double b, const double h )
	{
		Rinv = Eigen::MatrixXd::Zero( channels, channels );
		
        double hh12 = h * h / 12.0;
			
		xvec[0] = a;
		xvec[NPoints - 1] = b;
		detRvec[0] = std::numeric_limits<double>::max();
		detRvec[NPoints - 1] = std::numeric_limits<double>::max();
		
		double x = b - h;
		for ( int i = NPoints - 1; i > 0; i--, x -= h )
		{
			fill_V( x );

			Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
			U = 12.0 * Wmat.inverse() - 10.0 * I;
			R = U - Rinv;
			Rinv = R.inverse();

			xvec[i] =  x;
			detRvec[i] = R.determinant();	
		}	
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

		const double eps_abs = 1.0e-8;
		const double eps_rel = 1.0e-8;

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
	
	double adiabatic_potential( const double r, const int k )
	{
		fill_V( r );
		es.compute( V );
		
		return es.eigenvalues()(k); 
	}

	std::pair<double, double> find_classical_turning_points( const double Energy, double x_lb = 6.0, double x_rb = 22.0 )
	// если задать Energy = -1.0e-13 для этого потенциала, то поворотные точки получаются (6.0, 22.0), поэтому эти значения 
	// можно задать в качестве дефолтных, т.к. для любого связанного состояния поворотные точки будут лежать в этом интервале. 
	{
		double dx = 0.5;

		Eigen::VectorXd Veig_prev = Eigen::VectorXd::Zero(channels);

		std::function<double(double)> tmp;

		double tp;
		// left and right turning points for different channels
		std::vector<double> left_tps(channels);
		std::vector<double> right_tps(channels); 

		for ( double x = x_lb; x < x_rb; x += dx, Veig_prev = Veig )
		{
			fill_V( x );
			es.compute( V );
			Veig = es.eigenvalues();
			
			for ( int k = 0; k < channels; ++k )
			{
				if ( (Veig(k) - Energy) * (Veig_prev(k) - Energy) < 0 )
				{
					//std::cout << "Localized zero of " << k << " channel. Interval: " << x - dx << ", " << x << std::endl;
				
					if ( Veig(k) - Energy < 0 )
					{
						//std::cout << "Left turning point. Calling brent " << std::endl;
						tmp = [=]( const double r ) { return adiabatic_potential( r, k ) - Energy; }; 
										
						tp = brent( tmp, x - dx, x );
						//std::cout << "Turning point: " << tp << std::endl;

						left_tps[k] = tp;	
					}

					else if ( Veig(k) - Energy > 0 )
					{
						//std::cout << "Right turning point. Calling brent " << std::endl;
						tmp = [=]( const double r ) { return adiabatic_potential( r, k ) - Energy; };

						tp = brent( tmp, x - dx, x );
						//std::cout << "Turning point: " << tp << std::endl;

						right_tps[k] = tp;	
					}
				}
			}
		}

		double left_tp = *std::min_element( left_tps.begin(), left_tps.end() );
	   	double right_tp = *std::max_element( right_tps.begin(), right_tps.end() );

		return std::make_pair( left_tp, right_tp );
	}

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

	double D( const double a, const double b, const double h, const int i_match, const double E )
	{
		Eigen::MatrixXd Rm, Rmp1;

		propagateForward( E, a, h, i_match, Rm );
		propagateBackward( E, b, h, i_match, Rmp1 );
		
		return (Rm - Rmp1.inverse()).determinant();
	}

	void show_vector( std::string const & name, std::vector<double> const & v )
	{
		std::cout << name << ": ";
		for ( size_t k = 0; k < v.size(); ++k )
			std::cout << v[k] << " ";
		std::cout << "; size: " << v.size() << std::endl;
	}

    void show_vector( std::string const & name, std::vector<std::pair<double, double>> const & v )
    {
        std::cout << name << ": ";
        for ( size_t k = 0; k < v.size(); ++k )
            std::cout << "(" << v[k].first << ", " << v[k].second << ") ";
        std::cout << "; size: " << v.size() << std::endl;
    }

    void show_vector( std::string const & name, std::vector<Interval> const & v )
    {
        std::cout << "----" << std::endl;
        std::cout << name << ": " << std::endl;
        for ( size_t k = 0; k < v.size(); ++k )
            std::cout << "interval(E_min = " << v[k].min << ", E_max = " << v[k].max << ", nodes_count_min = " << v[k].node_count_min << ", nodes_count_max = " << v[k].node_count_max << ")" << std::endl;
        std::cout << "---- size: " << v.size() << std::endl;
    }

	int find_i_match( const double a, const double b, const double h, std::vector<double> const & nodesPos )
	// кладем точку матчинга в середину наибольшего интервала между узлами и границами отрезка; 
	// ищем наибольший интервал между узлами, координаты которых лежат в векторе nodesPos
	// также этот интервал может быть между началом отрезка, a, и первым узлом 
	// или между последним узлом и концом отрезка b.
	{
		double pos; // выбранная координата для точки матчинга
		double delta, new_delta; // длина отрезка, в середине которого будет лежать точка матчинга

		// создадим вспомогательный вектор, добавим к списку узлов концы отрезка
		std::vector<double> tmp( nodesPos.size() + 2 );
		std::copy( nodesPos.begin(), nodesPos.end(), tmp.begin() + 1 );
		tmp[0] = a;
		tmp.back() = b;

		//show_vector( "nodesPos", nodesPos );
		//show_vector( "tmp", tmp );

		// положим в начале точку матчинга в интервал между первой парой точек
		pos = 0.5 * (tmp[0] + tmp[1]);
		delta = tmp[1] - tmp[0];

		// затем проходим по следующим парам точек, 
		// если расстояние между ними больше чем текущее наибольшее, то запоминаем его как новое наибольшее
		size_t size_ = tmp.size();	
		for ( size_t k = 2; k < size_; ++k )
		{
			new_delta = tmp[k] - tmp[k - 1];

			if ( new_delta > delta )
			{
				delta = new_delta;
				pos = tmp[k - 1] + 0.5 * delta;	
			}
		}

		// находим координату точки матчинга
		int i_match = NPoints;
		for ( double approx = b; approx > pos; approx -= h, --i_match );

		return i_match;
	}

    Eigenvalue locate_eigenvalue_detR( double E_min, double E_max, const int nodes1, const int nodes2, const double a, const double b, const double h, const double eps )
    // nodes1 -- количество узлов при энергии E_min
    // nodes2 -- количество узлов при энергии E_max
    // предполагаем, что nodes2 = nodes + 1
    // eps: погрешность нахождения собственнного значения
    {
			//std::cout << "(locate_eigenvalue_detR) is called with parameters: E_min = " << E_min << "; E_max = " << E_max << "; nodes1 = " << nodes1 << "; nodes2 = " << nodes2 << std::endl;

        std::vector<double> nodes;
            
        int nodes_current = 0;
        double E_mean = 0;

        while ( std::abs(E_max - E_min) > eps )
        {
            E_mean = 0.5 * (E_max + E_min);

            propagateDetR( E_mean, a, b, h, nodes );
            nodes_current = nodes.size();
			//show_vector( "(locate_eigenvalue_detR) nodes", nodes );
			//std::cout << "(locate_eigenvalue_detR) E_mean: " << E_mean << "; nodes_current: " << nodes_current << std::endl;

			if ( nodes_current < nodes1 )
			{
				std::cerr << "(locate_eigenvalue_detR) Wrong behaviour of detR function! Number of nodes decreased!" << std::endl;
				exit( 1 );
			}

            if ( nodes_current == nodes1 )
            {
                E_min = E_mean;
				//std::cout << "Nodes current: " << nodes_current << "; nodes1: " << nodes1 << "; nodes2: " << nodes2 << "; updating E_min = " << E_min << "; E_max = " << E_max << std::endl;
            } 
            else if ( nodes_current == nodes2 )
            {
                E_max = E_mean;
				//std::cout << "Nodes current: " << nodes_current << "; nodes1: " << nodes1 << "; nodes2: " << nodes2 << "; E_min = " << E_min << "; updating E_max to " << E_max << std::endl;
            }
        }
        
		Eigenvalue eig( E_min, E_max, nodes1, nodes2 );
        return eig;
    }

    std::vector<Interval> separate_eigenvalues( double E_min, double E_max, const int nodes_min, const int nodes_max, const double a, const double b, const double h, const double eps )
    // nodes_min -- количество узлов функции |R|(x) при E = E_min
    // nodes_max -- количество узлов функции |R|(x) при E = E_max
    {
        std::cout << "(separate_eigenvalues) is called with parameters E_min: " << E_min << "; E_max: " << E_max << "; nodes_min: " << nodes_min << "; nodes_max: " << nodes_max << std::endl;
        std::vector<Interval> intervals;
        std::vector<Interval> recursive_intervals;

        std::vector<double> nodes;
        int nodes_current = 0;
        double E_mean = 0;
    
        if ( std::abs(E_min - E_max) > eps )
        {
            E_mean = 0.5 * (E_max + E_min);
            propagateDetR( E_mean, a, b, h, nodes );
            nodes_current = nodes.size();

            std::cout << "(separate_eigenvalues) E_mean: " << E_mean << "; nodes_current: " << nodes_current << std::endl; 

            // если отделили один узел с левого конца интервала
            if ( nodes_current == nodes_min + 1 )
            {
                intervals.emplace_back( E_min, E_mean, nodes_min, nodes_current );
                std::cout << "(separate_eigenvalues) isolated eigenvalue from the beginning of the interval! pushing interval " << E_min << ", " << E_mean << std::endl;
                
                if ( nodes_max - nodes_min == 2 )
                {
                    intervals.emplace_back( E_mean, E_max, nodes_current, nodes_max );
                    std::cout << "(separate_eigenvalues) pushing interval: " << E_mean << ", " << E_max << std::endl;
                    std::cout << "(separate_eigenvalues) all eigenvalues are separated!" << std::endl;
                }
                else
                {
                    std::cout << "(separate_eigenvalues) making recursive call!" << std::endl;
                    recursive_intervals = separate_eigenvalues( E_mean, E_max, nodes_current, nodes_max, a, b, h, eps );
                    std::cout << "(separate_eigenvalues) returned recursive_intervals: " << std::endl;
                    show_vector( "recursive_intervals", recursive_intervals );
                    intervals.insert( intervals.end(), recursive_intervals.begin(), recursive_intervals.end() );
                }
            }
            // если отделили некоторое количество узлов с левого конца
            else if ( nodes_current != nodes_min ) 
            {
                std::cout << "(separate_eigenvalues) making recursive call!" << std::endl;
                recursive_intervals = separate_eigenvalues( E_min, E_mean, nodes_min, nodes_current, a, b, h, eps );
                std::cout << "(separate_eigenvalues) returned recursive_intervals: " << std::endl;
                show_vector( "recursive_intervals", recursive_intervals );
                intervals.insert( intervals.end(), recursive_intervals.begin(), recursive_intervals.end() );

                if ( nodes_current != nodes_max )
                {
                    std::cout << "(separate eigenvalues) making recursive call!" << std::endl;
                    recursive_intervals = separate_eigenvalues( E_mean, E_max, nodes_current, nodes_max, a, b, h, eps );
                    std::cout << "(separate eigenvalues) return recursive intervals: " << std::endl;
                    show_vector( "recursive_intervals", recursive_intervals );
                    intervals.insert( intervals.end(), recursive_intervals.begin(), recursive_intervals.end() );
                }
            }
            // если таким образом не отделили узлов с левого конца
            else 
            {
                std::cout << "(separate_eigenvalues) making recursive call!" << std::endl;
                recursive_intervals = separate_eigenvalues( E_mean, E_max, nodes_current, nodes_max, a, b, h, eps );
                std::cout << "(separate_eigenvalues) returned recursive_intervals: " << std::endl;
                show_vector( "recursive_intervals", recursive_intervals );
                intervals.insert( intervals.end(), recursive_intervals.begin(), recursive_intervals.end() );
            }
        }
        else
        {
            std::cout << "Degenerate eigenvalue!" << std::endl;
        }

        return intervals;
    }
        
    std::vector<Eigenvalue> findEigenvalues_detR( const double E_min, const double E_max, const int reduce_step )
    {
        std::cout << std::fixed << std::setprecision(13);
        std::vector<Eigenvalue> eigs;

		double E_step = 1.0e-6;
		//double E_step = 3.0e-6;
		//double E_step = 7.0e-7;

		int nodes_prev = 0;
        int nodes_current = 0;

		// нужно пересчитывать a, b раз в некоторое время!
		// пока что просто зададим их какими-то большими, чтобы точно хватило
		double a = 5.0;
		double b = 16.0;
		double h = (b - a) / (NPoints - 1);

		std::vector<double> nodesPos;
        const double eps = 1.0e-12;

        std::vector<Interval> intervals;
		Eigenvalue eig;

        // первоначальное обнаружение количество узлов при энергии E_min
        propagateDetR( E_min, a, b, h, nodesPos );
        nodes_prev = nodesPos.size();
        std::cout << "(findEigenvalues_detR) E_min: " << E_min << "; nodes (nodes_prev): " << nodes_prev << std::endl;

		int counter = 0;
		for ( double E = E_min; E <= E_max; E += E_step, nodes_prev = nodes_current, ++counter )
		{
			propagateDetR( E, a, b, h, nodesPos );
			nodes_current = nodesPos.size(); 

			std::cout << "(findEigenvalues_detR) E: " << E << "; nodes_current: " << nodes_current << std::endl;
           
			if ( nodes_current < nodes_prev )
			{
				std::cerr << "(findEigenvalues) Wrong behaviour of detR function! Number of nodes decreased!" << std::endl;
				exit( 1 );
			}

			if ( nodes_current - nodes_prev == 1 )
			{
				std::cout << "(findEigenvalues_detR) Non-degenerate eigenvalue is detected." << std::endl;
	            eig = locate_eigenvalue_detR( E - E_step, E, nodes_prev, nodes_current, a, b, h, eps );
                std::cout << "(findEigenvalues_detR) located eigenvalue: " << std::endl << eig << std::endl;

                eigs.push_back( eig );
            }
            else if ( nodes_current - nodes_prev > 1 )
            {
                intervals = separate_eigenvalues( E - E_step, E, nodes_prev, nodes_current, a, b, h, eps );
            
                std::cout << "(findEigenvalues_detR) intervals isolated." << std::endl;
				show_vector( "intervals", intervals ); 

                for ( size_t k = 0; k < intervals.size(); ++k )
                {
                    eig = locate_eigenvalue_detR( intervals[k].min, intervals[k].max, intervals[k].node_count_min, intervals[k].node_count_max, a, b, h, eps );
                    std::cout << "(findEigenvalues_detR) isolated eigenvalue: " << std::endl << eig << std::endl;

                    eigs.push_back( eig );
                }
            }
			
			if ( (counter % reduce_step == 0) && (counter != 0) )
			{
				E_step /= 2.0;
				std::cout << "(findEigenvalues_detR) step = " << counter << "; reducing E_step to " << E_step << std::endl; 
			}
        }

        return eigs;
    }

	std::vector<double> findEigenvalues( const double E_min, const double E_max )
	{
		std::cout << std::fixed << std::setprecision(13);
		std::vector<double> res;

		double E_step = 1.0e-5;

		size_t nodes_prev = 0;
		size_t nodes_current = 0;

		// нужно пересчитывать a, b раз в некоторое время!
		// пока что просто зададим их какими-то большими, чтобы точно хватило
		double a = 6.0;
		double b = 14.0;
		double h = (b - a) / (NPoints - 1);

		std::vector<double> nodesPos;
		int i_match = NPoints / 2; 

		std::function<double(double)> _D; // D(E), корень которой является собственным числом
		//double eps = 1.0e-6;

		for ( double E = E_min; E < E_max; E += E_step, nodes_prev = nodes_current )
		{
			propagateDetR( E, a, b, h, nodesPos );
			nodes_current = nodesPos.size(); 

			i_match  = find_i_match( a, b, h, nodesPos );

			std::cout << "E: " << E << "; nodes_current: " << nodes_current << "; i_match: " << i_match << std::endl;

			// а что если скакнуло больше чем на один?!!!
			if ( nodes_current - nodes_prev == 1 )
			{
				std::cout << "Non-degenerate eigenvalue is detected." << std::endl;
				std::cout << "Calling brent with: " << E - E_step << ", " << E << std::endl; 

				_D = [=]( const double E ) { return D(a, b, h, i_match, E); };			

				double eigvalue = brent( _D, E - E_step, E );
				std::cout << "eigvalue: " << eigvalue << std::endl;

				res.push_back( eigvalue );
			}
			//else if ( (nodes_current - nodes_prev > 1) && (nodes_prev > 0) ) // BE CAREFUL IF THE FIRST EIGVALUE IS DEGENERATE
			//{
				//std::cout << "Degenerate eigenvalue is probably detected." << std::endl;
				//std::cout << "Energy region: " << E - E_step << ", " << E << std::endl;		
				//std::cout << "Nodes current: " << nodes_current << "; nodes_prev: " << nodes_prev << std::endl;
				//std::cout << "Trying to locate eigenvalue by node counting procedure." << std::endl;

				//int degeneracy;
				//double eigvalue;
				//bool status = locateNodeCount( E - E_step, E, eigvalue, degeneracy, eps );

				//// PROCEDURE IS NOT FINISHED YET!
				//if ( status )
				//{
					//std::cout << "Degenerate eigenvalue is detected!" << std::endl;
					//std::cout << "eigvalue: " << eigvalue << "; degeneracy: " << degeneracy << std::endl;
				//}
				//// else ??
			//}

		}

		return res;
	}	

	double __D__( const double E )
	{	
		double a = 6.0;
		double b = 14.0;
		double h = (b - a) / (NPoints - 1);
		int i_match = NPoints / 2.0;
		//std::cout << "pos: " << a + i_match * h << std::endl;

		return D( a, b, h, i_match, E );
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
    int NPoints = 50000;

	int M = 0;
	int J = 3;

	Solver solver( channels, NPoints );
	solver.setAngularMomentum( J, M );

	//double Energy = -1.0e-10;
	//std::pair<double, double> tp = solver.find_classical_turning_points( Energy, 6.0, 22.0 );
	//std::cout << "TP: " << tp.first << " " << tp.second << std::endl;

    //int reduce_step = 10;
	//double E_min = -5.3125e-5; // 2.0e-3
	//double E_max = -5.0e-7;
	//double E_min = -6.625e-6;
	//double E_max = -5.0e-7;
	//double E_min = -5.5e-4;
    //double E_max = -5.39e-4;
    //solver.findEigenvalues_detR( E_min, E_max, reduce_step );	


    double Energy = -1.0e-10;
    double a = 5.0;
    double b = 45.0;
    double h = (b - a) / (NPoints - 1);
    std::vector<double> nodesPos;
    solver.propagateDetR( Energy, a, b, h, nodesPos ); 

    std::cout << "nodesPos.size(): " << nodesPos.size() << std::endl;

	// ----------------------------------------------------------------
	//double E_min = -5.5e-4;
	//double E_max = -5.3e-4;
	//double dE = 1.0e-6;

	//for ( double E = E_min; E < E_max; E += dE )
		//std::cout << E << " " << solver.__D__( E ) << std::endl;
	// ----------------------------------------------------------------

	// ----------------------------------------------------------------
	//double E_min = -5.5e-4;
	//double E_max = -5.3e-4;
	
	//solver.findEigenvalues( E_min, E_max );
	// ----------------------------------------------------------------

	// ----------------------------------------------------------------
	//double a = 6.0;
	//double b = 14.0;
	//double h = (b - a) / (NPoints - 1);
	//std::vector<double> xvec(NPoints);
	//std::vector<double> detRvec(NPoints);

	//double Energy = -0.0000505;

	//solver.propagateDetR_testing( Energy, xvec, detRvec, a, b, h );	

	//int nodes = 0;
	//for ( size_t k = 0; k < xvec.size(); ++k )
	//{
		//if ( detRvec[k] < 0.0 )
		//{
			//std::cout << xvec[k] << " " << detRvec[k] << std::endl;
			//++nodes;
		//}
	//}
	//std::cout << "Nodes: " << nodes << std::endl;

	//std::ofstream file1( "file1.txt" );
	//for ( size_t k = 0; k < xvec.size(); ++k )
			//file1 << xvec[k] << " " << detRvec[k] << std::endl;	

	//solver.propagateDetR_backward_testing( Energy, xvec, detRvec, a, b, h );

	//std::ofstream file2( "file2.txt" );
	//for ( size_t k = 0; k < xvec.size(); ++k )
		//file2 << xvec[k] << " " << detRvec[k] << std::endl;	
	// ----------------------------------------------------------------

	return 0;
}
