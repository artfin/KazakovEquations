#include <Eigen/Dense>

#include <iostream>
#include <iomanip>

int main()
{
	std::cout << std::fixed << std::setprecision(14);	

	Eigen::MatrixXd m = Eigen::MatrixXd::Zero(4, 4);
	m << -0.00205704183097, 0.00000000000000, 0.04641332685576, 0.00000000000000,
	      0.00000000000000, 0.04037738825849, 0.00000000000000, 0.04002397211674,
	 	  0.04641332685576, 0.00000000000000, 0.04272444356415, 0.00000000000000,
 		  0.00000000000000, 0.04002397211674, 0.00000000000000, 0.03967364830178;

	std::cout << "m: " << std::endl << m << std::endl;
	
	Eigen::FullPivLU<Eigen::MatrixXd> lu( m );
	lu.setThreshold( 1.0e-5 );

	std::cout << "Kernel: " << std::endl << lu.kernel() << std::endl;

	return 0;
}
