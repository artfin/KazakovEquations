#include "./preliminary_eigenvalue_finder.hpp"
    
void PreliminaryEigenvalueFinder::show_vector( std::string const & name, std::vector<Eigenvalue> const & v )
{
    std::cout << "----" << std::endl;
    std::cout << name << ": " << std::endl;
    for ( size_t k = 0; k < v.size(); ++k )
        std::cout << "interval(E_min = " << v[k].get_min() << ", E_max = " << v[k].get_max() << ", nodes_count_min = " << v[k].get_node_count_min() << ", nodes_count_max = " << v[k].get_node_count_max() << ")" << std::endl;
    std::cout << "---- size: " << v.size() << std::endl;
}

Eigenvalue PreliminaryEigenvalueFinder::convergeToEigenvalue( Eigenvalue const & e, const double eps )
// function assumes that eigenvalue is already isolated,
// meaning that node_count_max = node_count_min + 1
// function returns reduced energy interval for eigenvalue,
// that complies with epsilon value
{
    //std::cout << "(convergeToEigenvalue) e.min: " << e.min << "; e.max: " << e.max << std::endl;

    int nodes;
    double Energy = 0;

    double min = e.get_min();
    double max = e.get_max();

    int nc_min = e.get_node_count_min();
    int nc_max = e.get_node_count_max();

    double a, b, h;
    std::pair<double, double> tp = equations->interpolate( max, energy_dict );
    equations->calculate_boundaries( tp, &a, &b, &h );

    // отрезок очень небольшой, нет смысла пересчитывать пределы
    while ( max - min > eps )
    {
        Energy = 0.5 * (min + max);
        nodes = equations->countEigenvalues( parity, Energy, a, b, h );

        if ( nodes == nc_min )
            min = Energy;
        else if ( nodes == nc_max )
            max = Energy;
    }

    // copy the previous eigenvalue; update min and max fields
    Eigenvalue eig = e;
    eig.set_min( min );
    eig.set_max( max );
    eig.update_value();

    return eig;
}

std::vector<Eigenvalue> PreliminaryEigenvalueFinder::findEigenvalues( const double E_min, const double E_max )
{
    std::pair<double, double> tp = equations->interpolate( E_min, energy_dict );
    
    double a, b, h;
    equations->calculate_boundaries( tp, &a, &b, &h );
    
    int nodes_min = equations->countEigenvalues( parity, E_min, a, b, h );

    tp = equations->interpolate( E_max, energy_dict );
    equations->calculate_boundaries( tp, &a, &b, &h );
    
    int nodes_max = equations->countEigenvalues( parity, E_max, a, b, h );

    std::vector<Eigenvalue> res = recursive_find( E_min, E_max, nodes_min, nodes_max );

    for ( size_t k = 0; k < res.size(); ++k )
        res[k].setAngularMomentum( equations->get_J(), equations->get_M() );

    return res;
}

std::vector<Eigenvalue> PreliminaryEigenvalueFinder::recursive_find( const double E_min, const double E_max, const int nodes_min, const int nodes_max )
{
    std::vector<Eigenvalue> res, tmp;
    
    if ( nodes_max - nodes_min == 1 )
    {
        std::cout << "(findEigenvalues) a node is already isolated." << std::endl;
        res.emplace_back( parity, E_min, E_max, nodes_min, nodes_max );
        return res; 
    }

    std::cout << "(findEigenvalues) is called with E_min: " << E_min << "; E_max: " << E_max << "; nodes_min: " << nodes_min << "; nodes_max: " << nodes_max << std::endl;

    int found_eigenvalues = 0;
    int eigenvalues_to_find = nodes_max - nodes_min;
    //std::cout << "(findEigenvalues) eigenvalues to find: " << eigenvalues_to_find << std::endl;

    double E_mean;
    int nodes;


    std::pair<double, double> tp; // pair of turning points
    double a, b, h;

    while ( found_eigenvalues != eigenvalues_to_find )
    {
        E_mean = 0.5 * (E_min + E_max);
        tp = equations->interpolate( E_mean, energy_dict );  
        equations->calculate_boundaries( tp, &a, &b, &h );

        nodes = equations->countEigenvalues( parity, E_mean, a, b, h );

        std::cout << "(recursive_find) E_mean: " << E_mean << "; a: " << a << "; b: " << b << "; nodes: " << nodes << std::endl;

        // если отделили 1 узел с левого края
        if ( nodes == nodes_min + 1 )
        {
            res.emplace_back( parity, E_min, E_mean, nodes_min, nodes );
            //std::cout << "(findEigenvalues) 1 eigenvalue isolated on (" <<  E_min << ", " << E_mean << ")" << std::endl;
            ++found_eigenvalues;

            if ( eigenvalues_to_find == 2 )
            {
                res.emplace_back( parity, E_mean, E_max, nodes, nodes_max );
                //std::cout << "(findEigenvalues) 1 eigenvalue isolated on (" << E_mean << ", " << E_max << ")" << std::endl;
                ++found_eigenvalues;
            }
            else
            {
                tmp = recursive_find( E_mean, E_max, nodes, nodes_max );
                res.insert( res.end(), tmp.begin(), tmp.end() );

                //std::cout << "(findEigenvalues) eigenvalues added: " << tmp.size() << std::endl;
                found_eigenvalues += tmp.size();
            }
        } 
        // если отделелили некоторое количество узлов с левого края
        else if ( nodes - nodes_min > 1 )
        {
            tmp = recursive_find( E_min, E_mean, nodes_min, nodes );
            res.insert( res.end(), tmp.begin(), tmp.end() );

            //std::cout << "(findEigenvalues) eigenvalues added: " << tmp.size() << std::endl;
            found_eigenvalues += tmp.size();
            
            if ( nodes != nodes_max )
            {
                tmp = recursive_find( E_mean, E_max, nodes, nodes_max );
                res.insert( res.end(), tmp.begin(), tmp.end() );

                //std::cout << "(findEigenvalues) eigenvalues added: " << tmp.size() << std::endl;
                found_eigenvalues += tmp.size();
            }
        }
        // если не отделили ни одного узла с левого края    
        else
        {
            tmp = recursive_find( E_mean, E_max, nodes, nodes_max );
            res.insert( res.end(), tmp.begin(), tmp.end() );

            //std::cout << "(findEigenvalues) eigenvalues added: " << tmp.size() << std::endl;
            found_eigenvalues += tmp.size();
        }
    }  

    return res; 
}



