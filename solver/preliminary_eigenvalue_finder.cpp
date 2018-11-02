#include "./preliminary_eigenvalue_finder.hpp"
    
void PreliminaryEigenvalueFinder::show_vector( std::string const & name, std::vector<Eigenvalue> const & v )
{
    std::cout << "----" << std::endl;
    std::cout << name << ": " << std::endl;
    for ( size_t k = 0; k < v.size(); ++k )
        std::cout << "interval(E_min = " << v[k].min << ", E_max = " << v[k].max << ", nodes_count_min = " << v[k].node_count_min << ", nodes_count_max = " << v[k].node_count_max << ")" << std::endl;
    std::cout << "---- size: " << v.size() << std::endl;
}

Eigenvalue PreliminaryEigenvalueFinder::convergeToEigenvalue( Eigenvalue const & e, const double a, const double b, const double h, const double eps )
// function assumes that eigenvalue is already isolated,
// meaning that node_count_max = node_count_min + 1
// function returns reduced energy interval for eigenvalue,
// that complies with epsilon value
{
    //std::cout << "(convergeToEigenvalue) e.min: " << e.min << "; e.max: " << e.max << std::endl;

    int nodes;
    double Energy = 0;

    double min = e.min;
    double max = e.max;

    int nc_min = e.node_count_min;
    int nc_max = e.node_count_max;

    while ( max - min > eps )
    {
        Energy = 0.5 * (min + max);
        nodes = equations->countEigenvalues( Energy, a, b, h );

        if ( nodes == nc_min )
            min = Energy;
        else
            max = Energy;
    }

    return Eigenvalue( min, max, nc_min, nc_max );
}

std::vector<Eigenvalue> PreliminaryEigenvalueFinder::findEigenvalues( const double E_min, const double E_max, const double a, const double b, const double h )
{
    int nodes_min = equations->countEigenvalues( E_min, a, b, h );
    int nodes_max = equations->countEigenvalues( E_max, a, b, h );

    return recursive_find( E_min, E_max, nodes_min, nodes_max, a, b, h );
}

std::vector<Eigenvalue> PreliminaryEigenvalueFinder::recursive_find( const double E_min, const double E_max, const int nodes_min, const int nodes_max, const double a, const double b, const double h )
{
    //std::cout << "(findEigenvalues) is called with E_min: " << E_min << "; E_max: " << E_max << "; nodes_min: " << nodes_min << "; nodes_max: " << nodes_max << std::endl;

    int found_eigenvalues = 0;
    int eigenvalues_to_find = nodes_max - nodes_min;
    //std::cout << "(findEigenvalues) eigenvalues to find: " << eigenvalues_to_find << std::endl;

    double E_mean;
    int nodes;

    std::vector<Eigenvalue> res, tmp;

    while ( found_eigenvalues != eigenvalues_to_find )
    {
        E_mean = 0.5 * (E_min + E_max);
        nodes = equations->countEigenvalues( E_mean, a, b, h );
   
        // если отделили 1 узел с левого края
        if ( nodes == nodes_min + 1 )
        {
            res.emplace_back( E_min, E_mean, nodes_min, nodes );
            //std::cout << "(findEigenvalues) 1 eigenvalue isolated on (" <<  E_min << ", " << E_mean << ")" << std::endl;
            ++found_eigenvalues;

            if ( eigenvalues_to_find == 2 )
            {
                res.emplace_back( E_mean, E_max, nodes, nodes_max );
                //std::cout << "(findEigenvalues) 1 eigenvalue isolated on (" << E_mean << ", " << E_max << ")" << std::endl;
                ++found_eigenvalues;
            }
            else
            {
                tmp = recursive_find( E_mean, E_max, nodes, nodes_max, a, b, h );
                res.insert( res.end(), tmp.begin(), tmp.end() );

                //std::cout << "(findEigenvalues) eigenvalues added: " << tmp.size() << std::endl;
                found_eigenvalues += tmp.size();
            }
        } 
        // если отделелили некоторое количество узлов с левого края
        else if ( nodes - nodes_min > 1 )
        {
            tmp = recursive_find( E_min, E_mean, nodes_min, nodes, a, b, h );
            res.insert( res.end(), tmp.begin(), tmp.end() );

            //std::cout << "(findEigenvalues) eigenvalues added: " << tmp.size() << std::endl;
            found_eigenvalues += tmp.size();
            
            if ( nodes != nodes_max )
            {
                tmp = recursive_find( E_mean, E_max, nodes, nodes_max, a, b, h );
                res.insert( res.end(), tmp.begin(), tmp.end() );

                //std::cout << "(findEigenvalues) eigenvalues added: " << tmp.size() << std::endl;
                found_eigenvalues += tmp.size();
            }
        }
        // если не отделили ни одного узла с левого края    
        else
        {
            tmp = recursive_find( E_mean, E_max, nodes, nodes_max, a, b, h );
            res.insert( res.end(), tmp.begin(), tmp.end() );

            //std::cout << "(findEigenvalues) eigenvalues added: " << tmp.size() << std::endl;
            found_eigenvalues += tmp.size();
        }
    }  

   return res; 
}
