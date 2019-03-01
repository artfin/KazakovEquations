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
            auto& fn = *static_cast<std::function<double(double)>*>(vf);
            return fn(d);
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

    equations->propagateForwardToMatch( E, a, h, i_match, Rm );
    equations->propagateBackwardToMatch( E, b, h, i_match, Rmp1 );

    //std::cout << "(D) i_match: " << i_match << "; Rm: " << std::endl << Rm << std::endl << "; Rmp1_inv: " << std::endl << Rmp1.inverse() << std::endl;

    std::cout << std::fixed << std::setprecision(16);
    std::cout << "(D) i_match: " << i_match << "; a: " << a << "; b: " << b << "; E: " << E*constants::HTOCM << "; determinant: " << (Rm - Rmp1.inverse()).determinant() << std::endl;

    Eigen::MatrixXd D = Rm - Rmp1.inverse();
    //std::cout << "(D) D: " << std::endl << D << std::endl;
    return D.determinant();
}

double PreciseEigenvalueFinder::precise_eigenvalue_calculation( const int i_match, const double E1, const double E2, double & a, double & b, double & h )
{
    std::pair<double, double> tp = equations->interpolate( E2, energy_dict );  
    equations->calculate_boundaries( tp, &a, &b, &h );
    // энергии E1, E2 уже достаточно близки, нет никакого смысла пересчитывать a, b, h для дальнейших шагов

    std::function<double(double)> __D__ = [=]( const double E ) { return D(a, b, h, i_match, E); };
    return brent( __D__, E1, E2 ); 
}

int PreciseEigenvalueFinder::effecient_i_match( const double E1, const double E2, const int NPoints, const int parts )
{
    int len = parts - 1; // only internal dots 

    // энергии очень близки, достаточно одного набора поворотных точек
    double a, b, h;
    std::pair<double, double> tp = equations->interpolate( E2, energy_dict );
    equations->calculate_boundaries( tp, &a, &b, &h );

    int step = NPoints / parts;

    // ENERGY = E1!
    equations->propagateForwardFull( E1, a, h, Rm_vector1, step );
    equations->propagateBackwardFull( E1, b, h, Rmp1_vector1, step );

    // ENERGY = E2!
    equations->propagateForwardFull( E2, a, h, Rm_vector2, step );
    equations->propagateBackwardFull( E2, b, h, Rmp1_vector2, step );

    double D1, D2;
    int i_match = step;
    for ( int k = 0; k < len; ++k, i_match += step )
    {
        D1 = (Rm_vector1[k] - Rmp1_vector1[k]).determinant();
        D2 = (Rm_vector2[k] - Rmp1_vector2[k]).determinant();

        //std::cout << "(effecient_i_match) i_match: " << i_match << "; D1: " << D1 << "; D2: " << D2 << std::endl;
        if ( D1 * D2 < 0.0 )
            return i_match;
    }

    // если подходящая точка сшивки не найдена, то возвращаем -1.
    return -1;
}

int PreciseEigenvalueFinder::find_i_match( const double E1, const double E2, const int NPoints, const int parts )
{
    int i_match;
    double D1, D2;

    std::pair<double, double> tp = equations->interpolate( E2, energy_dict );
    
    double a, b, h;
    equations->calculate_boundaries( tp, &a, &b, &h );

    for ( int k = 1; k < parts; ++k )
    {
        i_match = NPoints / parts * k;
		//std::cout << "i_match: " << i_match << std::endl;
        
        D1 = D(a, b, h, i_match, E1);
        D2 = D(a, b, h, i_match, E2);

        if ( D1 * D2 < 0.0 )
        {
            //std::cerr << "D1 * D2 < 0.0! i_match = " << i_match << " is appropriate!" << std::endl;
            return i_match;
        }
    }

    return -1;
}


void PreciseEigenvalueFinder::reserve_space_for_search( const int parts )
{
    int len = parts - 1; 
    
    Rm_vector1.resize( len );
    Rmp1_vector1.resize( len );
   
    Rm_vector2.resize( len );
    Rmp1_vector2.resize( len );
}

