#include "./equations.hpp"

Equations::Equations( int channels, int NPoints ) : channels(channels), NPoints(NPoints)
{
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

void Equations::fill_V( const double r )
{
	// заполняем элементы вектора W
	fill_W_elements( r );

	for ( int P = 0; P < V.rows(); ++P )
		for ( int Lprime = 0; Lprime < V.cols(); ++Lprime )
		//V(P, Lprime) = compute_W_sum( 2*P, 2*Lprime ) + delta( 2*P, 2*Lprime ) * (hbar*hbar*2.0*P*(2.0*P+1.0)/2.0/Inten + hbar*hbar/2.0/mu/r/r*(J*(J+1) + 2.0*P*(2.0*P+1.0)));
			V(P, Lprime) = compute_W_sum( P, Lprime ) + delta(P, Lprime) * (hbar*hbar*P*(P+1)/2.0/Inten + hbar*hbar/2.0/mu/r/r*(J*(J+1) + P*(P+1)));
}

void Equations::propagateDetR( const double Energy, const double a, const double b, const double h, std::vector<double> & nodesPos )
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

void Equations::propagateForward( const double Energy, const double a, const double h, const int i_match, Eigen::MatrixXd & resRm, bool save ) 
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

void Equations::propagateBackward( const double Energy, const double b, const double h, const int i_match, Eigen::MatrixXd & resRmp1, bool save ) 
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

