#include "./equations.hpp"

Equations::Equations( int channels_, int NPoints ) : channels(channels_ / 2), NPoints(NPoints) 
{
    std::cout << "(equations constructor) Effective number of channels: " << 2 * channels << "; real size of matrices: " << channels << std::endl;
    W.resize( 6 );

	V = Eigen::MatrixXd::Zero( channels, channels );
    I = Eigen::MatrixXd::Identity( channels, channels );
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

void Equations::setAngularMomentum( const int J, const int M )
{
	this->J = J;
	this->M = M;
}

void Equations::fill_W_elements( const double R )
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

double Equations::angularMatrixElements( const int P, const int l, const int Lprime )
// gsl_sf_coupling_3j принимает удвоенные значения спинов и проекций спинов
// возвращает значение матричного элемента <PM|P_l|L'M>
{
	return std::pow(-1.0, M) * std::sqrt( (2.0*P+1.0) * (2.0*Lprime+1.0) ) * gsl_sf_coupling_3j(2*Lprime, 2*l, 2*P, 0, 0, 0 ) * gsl_sf_coupling_3j(2*Lprime, 2*l, 2*P, 2*M, 0, -2*M );
}
	
double Equations::compute_W_sum( const int P, const int Lprime)
// возвращает значение суммы \sum_l W_l * <PM|P_l|L'M>
{
	double result = 0.0;
	for ( size_t i = 0; i < W.size(); ++i )
		result += W[i] * angularMatrixElements( P, 2 * i, Lprime );

	return result;
}

void Equations::fill_V( const double r, const Parity parity )
{
    assert( J >= 0 && "Please, set angular momentum using setAngularMomentum method!" );

	// заполняем элементы вектора W
	fill_W_elements( r );

    if ( parity == Parity::ODD )
    {
        for ( int P = 0; P < V.rows(); ++P )
        {
            for ( int Lprime = 0; Lprime < V.cols(); ++Lprime )
            {
                V(P, Lprime) = compute_W_sum(2*P+1, 2*Lprime+1) + delta(2*P+1, 2*Lprime+1) * (hbar*hbar*(2.0*P+1.0)*(2.0*P+2.0)/2.0/Inten + hbar*hbar/2.0/mu/r/r*(J*(J+1.0) + (2.0*P+1.0)*(2.0*P+2.0)));
            }
        }
    }
    else if ( parity == Parity::EVEN )
    {
        for ( int P = 0; P < V.rows(); ++P )
        {
            for ( int Lprime = 0; Lprime < V.cols(); ++Lprime )
            {
                V(P, Lprime) = compute_W_sum(2*P, 2*Lprime) + delta(2*P, 2*Lprime) * (hbar*hbar*2.0*P*(2.0*P+1.0)/2.0/Inten + hbar*hbar/2.0/mu/r/r*(J*(J+1.0) + 2.0*P*(2.0*P+1.0)));
            }
        }
    }
    else
    {
        for ( int P = 0; P < V.rows(); ++P )
        {
            for ( int Lprime = 0; Lprime < V.cols(); ++Lprime )
            {
                V(P, Lprime) = compute_W_sum( P, Lprime ); // + delta(P, Lprime) * (hbar*hbar*P*(P+1)/2.0/Inten + hbar*hbar/2.0/mu/r/r*(J*(J+1) + P*(P+1)));
            }
        }
    }
}

void Equations::propagateDetR( const Parity parity, const double Energy, const double a, const double b, const double h, std::vector<double> & nodesPos )
{   
    // clearing nodes vector!
    nodesPos.clear();

	Rinv = Eigen::MatrixXd::Zero( channels, channels );
		
    double hh12 = h * h / 12.0;
	double detR = std::numeric_limits<double>::max();
	double x = a + h;
	
    for ( int i = 1; i < NPoints - 1; i++, x += h )
	{
		fill_V( x, parity );

		Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
		U = 12.0 * Wmat.inverse() - 10.0 * I;
		R = U - Rinv;
		Rinv = R.inverse();

		detR = R.determinant();
		if ( detR < 0 )
			nodesPos.push_back( x );
	}
}

void Equations::propagateForwardFull( const Parity parity, const double Energy, const double a, const double h, std::vector<Eigen::MatrixXd> & Rm_vector, const int step )
{
    double hh12 = h * h / 12.0;
    double x = a + h;

    Rinv = Eigen::MatrixXd::Zero( channels, channels );

    int counter = 0;
    for ( int i = 1; i < NPoints - 1; i++, x += h )
    {
        fill_V( x, parity );

		Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
		U = 12.0 * Wmat.inverse() - 10.0 * I;
		R = U - Rinv;
		Rinv = R.inverse();
    
        if ( (i % step == 0) && (i != 0) )
        {
            //std::cout << "(propagateForwardFull) i: " << i << std::endl;
            Rm_vector[counter] = R;
            ++counter;
        } 
    } 
} 

void Equations::propagateBackwardFull( const Parity parity, const double Energy, const double b, const double h, std::vector<Eigen::MatrixXd> & Rmp1_vector, const int step )
{
    //std::cout << "(propagateBackwardFull) Energy: " << Energy << "; b: " << b << "; h: " << h << std::endl;
    
    double hh12 = h * h / 12.0;
    double x = b - h;

    Rinv = Eigen::MatrixXd::Zero( channels, channels );

    size_t size_ = Rmp1_vector.size() - 1;
    int counter = 0; 
    for ( int i = NPoints - 2; i > 0; i--, x -= h )
    {
        fill_V( x, parity );

        Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
        U = 12.0 * Wmat.inverse() - 10.0 * I;
        R = U - Rinv;
        Rinv = R.inverse();
        
        if ( ((i - 1) % step == 0) && (i - 1 != 0) ) 
        {
            //std::cout << "(propagateBackwardFull) i: " << i << std::endl;
            Rmp1_vector[size_ - counter] = Rinv;
            ++counter;
        }
    }	
}

void Equations::propagateForward( const Parity parity, const double Energy, const double a, const double h, const int i_match, Eigen::MatrixXd & resRm, bool save ) 
{
    double hh12 = h * h / 12.0;
    double x = a + h;

    Rinv = Eigen::MatrixXd::Zero( channels, channels );

    for ( int i = 1; i <= i_match; i++, x += h )
    {
        //std::cout << "(propagateForward) i: " << i << "; x: " << x << std::endl; 
        fill_V( x, parity );

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


void Equations::propagateBackward( Parity parity, const double Energy, const double b, const double h, const int i_match, Eigen::MatrixXd & resRmp1, bool save ) 
// resRmp1 = R_{m + 1}
{
    //std::cout << "(propagateBackward) Energy: " << Energy << "; b: " << b << "; h: " << h << std::endl;

    double hh12 = h * h / 12.0;
    double x = b - h;

    Rinv = Eigen::MatrixXd::Zero( channels, channels );

    for ( int i = NPoints - 2; i > i_match; i--, x -= h )
    {
        //std::cout << "(propagateBackward) i: " << i << "; x: " << x << std::endl;
        fill_V( x, parity );

        Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar; 
        U = 12.0 * Wmat.inverse() - 10.0 * I;
        R = U - Rinv;
        Rinv = R.inverse();

        if ( save )
        {
            Rinv_vector[i] = Rinv;
            Winv_vector[i] = Wmat.inverse();
        }
    }	

    resRmp1 = R;	
}

int Equations::countEigenvalues( const Parity parity, const double Energy, const double a, const double b, const double h )
{
    Rinv = Eigen::MatrixXd::Zero( channels, channels );

    double hh12 = h * h / 12.0;

    int n = 0;
    double x = a + h;
    for ( int i = 1; i < NPoints - 1; i++, x += h )
    {
        fill_V( x, parity );

        Wmat = I + hh12 * (Energy * I - V) * 2.0 * mu / hbar / hbar;
        U = 12.0 * Wmat.inverse() - 10.0 * I;
        R = U - Rinv;

        Rinv = R.inverse();

        n += countNegativeDiagonalElements();
    }

    return n;
} 

int Equations::countNegativeDiagonalElements()
{
    es.compute( R ); 
    eigenvalues = es.eigenvalues();

    int counter = 0;
    for ( int k = 0; k < channels; ++k )
        if ( eigenvalues(k) < 0.0 )
            ++counter;

    return counter;
}

double Equations::brent( std::function<double(double)> f, double xb1, double xb2, const double eps )
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

    do
    {
        iter++;

        status = gsl_root_fsolver_iterate( s );
        r = gsl_root_fsolver_root( s );

        x_lo = gsl_root_fsolver_x_lower( s );
        x_hi = gsl_root_fsolver_x_upper( s );
        status = gsl_root_test_interval( x_lo, x_hi, eps, eps );

        //if ( status == GSL_SUCCESS )
        //printf( "Converged:\n");

        //printf( "%5d [%.12f, %.12f] %.12f %.12f\n",
        //iter, x_lo, x_hi, r, x_hi - x_lo );	
    }
    while ( status == GSL_CONTINUE && iter < max_iter );	

    gsl_root_fsolver_free( s );

    return r;
}

Eigen::VectorXd Equations::adiabatic_potential_vector( const Parity parity, const double r )
{
    fill_V( r, parity );
    es.compute( V );
    return es.eigenvalues();
}

double Equations::adiabatic_potential_component( const Parity parity, const double r, const int k )
{
    fill_V( r, parity );
    es.compute( V );

    return es.eigenvalues()(k); 
}

std::pair<double, double> Equations::find_classical_turning_points( const Parity parity, const double Energy, double x_lb, double x_rb, const double eps )
{
    double dx = 0.5;

    Eigen::VectorXd Veig;
    Eigen::VectorXd Veig_prev = Eigen::VectorXd::Zero(channels);

    std::function<double(double)> tmp;

    double tp;
    // left and right turning points for different channels
    std::vector<double> left_tps(channels);
    std::vector<double> right_tps(channels); 

    for ( double x = x_lb; x < x_rb; x += dx, Veig_prev = Veig )
    {
        fill_V( x, parity );
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
                    tmp = [=]( const double r ) { return adiabatic_potential_component( parity, r, k ) - Energy; }; 

                    tp = brent( tmp, x - dx, x, eps );
                    //std::cout << "Turning point: " << tp << std::endl;

                    left_tps[k] = tp;
                }

                else if ( Veig(k) - Energy > 0 )
                {
                    //std::cout << "Right turning point. Calling brent " << std::endl;
                    tmp = [=]( const double r ) { return adiabatic_potential_component( parity, r, k ) - Energy; };

                    tp = brent( tmp, x - dx, x, eps );
                    //std::cout << "Turning point: " << tp << std::endl;

                    right_tps[k] = tp;
                }
            }
        }
    }

    // removing zero values for some channels
    left_tps.erase(
        std::remove( left_tps.begin(), left_tps.end(), 0.0 ),
        left_tps.end()
    );

    right_tps.erase(
        std::remove( right_tps.begin(), right_tps.end(), 0.0),
        right_tps.end()
    );


    double left_tp = *std::min_element( left_tps.begin(), left_tps.end() );
    double right_tp = *std::max_element( right_tps.begin(), right_tps.end() );

    return std::make_pair( left_tp, right_tp );
}

std::map<double, std::pair<double, double>> Equations::create_energy_dict( const Parity parity, const double E_min, const double E_max, const int energy_intervals, const double x_lb, const double x_rb, const double eps )
// x_lb, x_rb -- задают интервал, в котором будут искаться поворотные точки
// eps -- погрешность нахождения поворотной точки
{
    double lg_E_min = std::log10( -E_min );
    double lg_E_max = std::log10( -E_max );
    double energy_step = (lg_E_min - lg_E_max) / (energy_intervals - 1); 

    std::pair<double, double> tp;
    std::map<double, std::pair<double, double>> energy_dict;  

    std::vector<double> zero_tps;

    double lge = lg_E_max;
    for ( int counter = 0; counter < energy_intervals; lge += energy_step, counter++ )
    {
        double e = -std::pow(10.0, lge);
        tp = find_classical_turning_points( parity, e, x_lb, x_rb, eps );
        energy_dict[e] = tp;

        if ( (tp.first == 0) || (tp.second == 0) )
            zero_tps.push_back( e );
    } 
    
    //for ( auto it = energy_dict.begin(); it != energy_dict.end(); ++it )
    //std::cout << it->first << " " << it->second.first << " " << it->second.second << std::endl;

    // может быть такая ситуация, что при некоторой низкой энергии поворотных точек нет
    // в таком случае find_classical_turning_points(...) вернет пару (0.0, 0.0).
    // Однако в таком случае при какой-то последующей энергии у нас будут ненулевые поворотные точки
    // Предлагаю использовать их в в качестве поворотных и для этой энергии 
    // 
    // 1) Сначала надо найти первую снизу энергию, для которой поворотные точки ненулевые:
    double Elb_reference;
    std::pair<double, double> plb;
    for ( auto it = energy_dict.begin(); it != energy_dict.end(); ++it )
    {
        if ( (it->second.first != 0.0) && (it->second.second != 0.0) )
        {   
            plb = it->second;
            Elb_reference = it->first;
            break;
        }
    }
  
    // 2) Находим последнюю энергию, для которой поворотные точки ненулевые:
    double Eub_reference;
    std::pair<double, double> pub; 
    for ( auto it = energy_dict.rbegin(); it != energy_dict.rend(); ++it )
    {
        if ( (it->second.first != 0.0) && (it ->second.second != 0.0) )
        {
            pub = it->second;
            Eub_reference = it->first;
            break;
        }
    }

    // 3) Обходим те энергии, которые содержат нулевые пары и присваиваем им ненулевые 
    //    значения поворотных точек
    for ( size_t k = 0; k < zero_tps.size(); ++k )
    {
        if ( zero_tps[k] < Elb_reference )
        {
            auto it = energy_dict.find( zero_tps[k] );
            it->second = plb;
        }
        else if ( zero_tps[k] > Eub_reference )
        {
            auto it = energy_dict.find( zero_tps[k] );
            it->second = pub; 
        }
        else
        {
            std::cerr << "Problems with filling energy dict!" << std::endl;
            exit( 1 );
        }
    } 

    // 4) Проверяем, что правые поворотные точки идут в порядке возрастания
    auto it = energy_dict.begin();
    double prev = it->second.second;
    double next;

    it++;
    for ( ; it != energy_dict.end(); it++ )
    {
        next = it->second.second;

        if ( prev > next )
        {
            it->second.second = prev;

            std::cerr << "Problems with energy dict! Turning points are not in strict increasing order!" << std::endl;
            //exit( 1 );
        }
        else
        {
            prev = next;
        }
    }

    return energy_dict;
} 

std::pair<double, double> Equations::interpolate( const double e, const std::map<double, std::pair<double, double>> & energy_dict )
{
    double lower_bound = 0.0;
    double upper_bound = 0.0;
    const double eps = 1.0e-9;

    for(  auto it = energy_dict.begin(); it != energy_dict.end(); ++it )
    {
        if ( std::abs(it->first - e) < eps )
        {
            return it->second; 
        }

        if ( it->first > e )
        {
            upper_bound = it->first;
            break;
        }
    }

    for ( auto it = energy_dict.rbegin(); it != energy_dict.rend(); ++it )
    {
        if ( it->first < e )
        {
            lower_bound = it->first;
            break;
        }
    }

    if ( (lower_bound == 0.0) || (upper_bound == 0.0) )
    {
        std::cout << "Energy: " << e << std::endl;
        std::cerr << "Interpolation is not possible!" << std::endl;
        exit( 1 );
    }

    //std::cout << "(interpolate) lower_bound: " << lower_bound << "; upper_bound: " << upper_bound << std::endl;

    double lge = std::log10( -e );
    double lglb = std::log10( -lower_bound ); 
    double lgub = std::log10( -upper_bound );

    std::pair<double, double> tp1 = energy_dict.find(lower_bound)->second;
    std::pair<double, double> tp2 = energy_dict.find(upper_bound)->second;

    double tp1lb = tp1.first;
    double tp1ub = tp2.first;
    double tp2lb = tp1.second;
    double tp2ub = tp2.second;
       
    std::pair<double, double> res; 
    res.first = tp1lb + (lglb - lge) * (tp1lb - tp1ub) / (lgub - lglb);
    res.second = tp2lb + (lglb - lge) * (tp2lb - tp2ub) / (lgub - lglb); 

    return res;
} 

void Equations::calculate_boundaries( std::pair<double, double> const & tp, double * a, double * b, double * h )
{
    *a = tp.first - 1.5; 
    *b = tp.second + 6.0; 
    *h = ( *b - *a ) / (NPoints - 1);
   //std::cout << "(boundaries) a: " << *a << "; b: " << *b << std::endl;
} 

