#include "./precise_eigenvalue_finder.hpp"

double PreciseEigenvalueFinder::brent( std::function<double(double)> f, double xb1, double xb2 )
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

double PreciseEigenvalueFinder::D( const double a, const double b, const double h, const int i_match, const double E )
{
    Eigen::MatrixXd Rm, Rmp1;

    equations->propagateForward( E, a, h, i_match, Rm );
    equations->propagateBackward( E, b, h, i_match, Rmp1 );

    return (Rm - Rmp1.inverse()).determinant();
}

double PreciseEigenvalueFinder::precise_eigenvalue_calculation( const double a, const double b, const double h, const int i_match, const double E1, const double E2 )
{
    std::function<double(double)> __D__ = [=]( const double E ) { return D(a, b, h, i_match, E); };
    return brent( __D__, E1, E2 ); 
}

int PreciseEigenvalueFinder::find_i_match( const double a, const double b, const double h, const double E1, const double E2, const int parts )
{
    int i_match;
    double D1, D2;

    for ( int k = 1; k < parts; ++k )
    {
        i_match = NPoints / parts * k;
		//std::cout << "i_match: " << i_match << std::endl;
        
        D1 = D(a, b, h, i_match, E1);
        D2 = D(a, b, h, i_match, E2);

        if ( D1 * D2 < 0.0 )
        {
			std::cerr << "D1 * D2 < 0.0! i_match = " << i_match << " is appropriate!" << std::endl;
            return i_match;
        }
    }

    std::cerr << "i_match is not found!" << std::endl;
    exit( 1 );
}
