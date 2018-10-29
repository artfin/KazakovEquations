#include "./eigenfunction_finder.hpp"

void EigenfunctionFinder::calculate_eigenfunction( std::vector<Eigen::VectorXd> & eigfunc, const double E, const double a, const double b, const double h, int i_match  )
{
    std::cout << std::fixed << std::setprecision(14);

    Eigen::MatrixXd Rm, Rmp1;
    equations->propagateForward( E, a, h, i_match, Rm, true );
    equations->propagateBackward( E, b, h, i_match, Rmp1, true );
    Eigen::MatrixXd diff = Rm - Rmp1.inverse();

    //std::cout << "diff: " << std::endl << diff << std::endl << std::endl;

    Eigen::FullPivLU<Eigen::MatrixXd> lu( diff );
    lu = lu.setThreshold( 1.0e-5 );

    Eigen::MatrixXd ker = lu.kernel();
    //std::cout << "kernel: " << std::endl << ker << std::endl;

    Eigen::VectorXd fn_prev, fn, fn_next, psin;
    Eigen::VectorXd fn_matching_point;

    fn_matching_point = ker;
    fn_next = fn_matching_point;
    fn_prev = fn_matching_point;

    //std::cout << "---------------------------------------------" << std::endl << std::endl;
    //std::cout << "(calculate_eigenfunction)" << std::endl;

    for ( int k = i_match - 1;  k > 0; --k, fn_next = fn )
    {
        //std::cout << "k: " << k << "; Rinv: " << std::endl << equations->Rinv_vector[k] << std::endl;
        fn = equations->Rinv_vector[k] * fn_next;
        psin = equations->Winv_vector[k] * fn;
        eigfunc[k] = psin;	

        //std::cout << "eigfunc[" << k << "]: " << std::endl << eigfunc[k] << std::endl;	
    }	

    for ( int k = i_match + 1; k < NPoints - 1; ++k, fn_prev = fn )
    {
        //std::cout << "k: " << k << "; Rinv: " << std::endl << Rinv_vector[k] << std::endl;
        fn = equations->Rinv_vector[k] * fn_prev;
        psin = equations->Winv_vector[k] * fn;
        eigfunc[k] = psin;	

        //std::cout << "eigfunc[" << k << "]: " << std::endl << eigfunc[k] << std::endl;	
    }

    psin = equations->Winv_vector[i_match] * fn_matching_point;
    eigfunc[i_match] = psin;

    //std::cout << "eigfunc[" << i_match << "]: " << std::endl << eigfunc[i_match] << std::endl;
}

