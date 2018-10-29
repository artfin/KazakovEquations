#include "./preliminary_eigenvalue_finder.hpp"
    
void PreliminaryEigenvalueFinder::show_vector( std::string const & name, std::vector<Eigenvalue> const & v )
{
    std::cout << "----" << std::endl;
    std::cout << name << ": " << std::endl;
    for ( size_t k = 0; k < v.size(); ++k )
        std::cout << "interval(E_min = " << v[k].min << ", E_max = " << v[k].max << ", nodes_count_min = " << v[k].node_count_min << ", nodes_count_max = " << v[k].node_count_max << ")" << std::endl;
    std::cout << "---- size: " << v.size() << std::endl;
}

std::vector<Eigenvalue> PreliminaryEigenvalueFinder::separate_eigenvalues( double E_min, double E_max, const int nodes_min, const int nodes_max, const double a, const double b, const double h, const double eps )
// nodes_min -- количество узлов функции |R|(x) при E = E_min
// nodes_max -- количество узлов функции |R|(x) при E = E_max
{
    std::cout << "(separate_eigenvalues) is called with parameters E_min: " << E_min << "; E_max: " << E_max << "; nodes_min: " << nodes_min << "; nodes_max: " << nodes_max << std::endl;
    std::vector<Eigenvalue> intervals;
    std::vector<Eigenvalue> recursive_intervals;

    std::vector<double> nodes;
    int nodes_current = 0;
    double E_mean = 0;

    if ( std::abs(E_min - E_max) > eps )
    {
        E_mean = 0.5 * (E_max + E_min);
        equations->propagateDetR( E_mean, a, b, h, nodes );
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

Eigenvalue PreliminaryEigenvalueFinder::locate_eigenvalue( double E_min, double E_max, const int nodes1, const int nodes2, const double a, const double b, const double h, const double eps )
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

        equations->propagateDetR( E_mean, a, b, h, nodes );
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

std::vector<Eigenvalue> PreliminaryEigenvalueFinder::findEigenvalues( const double E_min, const double E_max, const double a, const double b, const double h, double E_step, const int reduce_step, const double eps )
{
	std::cout << std::fixed << std::setprecision(13);
	std::vector<Eigenvalue> eigs;

	int nodes_prev = 0;
	int nodes_current = 0;

	std::vector<double> nodesPos;

	std::vector<Eigenvalue> intervals;
	Eigenvalue eig;

	// первоначальное обнаружение количество узлов при энергии E_min
	equations->propagateDetR( E_min, a, b, h, nodesPos );
	nodes_prev = nodesPos.size();
	std::cout << "(findEigenvalues_detR) E_min: " << E_min << "; nodes (nodes_prev): " << nodes_prev << std::endl;

	int counter = 0;
	for ( double E = E_min; E < E_max; E += E_step, nodes_prev = nodes_current, ++counter )
	{
		equations->propagateDetR( E, a, b, h, nodesPos );
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
			eig = locate_eigenvalue( E - E_step, E, nodes_prev, nodes_current, a, b, h, eps );
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
				eig = locate_eigenvalue( intervals[k].min, intervals[k].max, intervals[k].node_count_min, intervals[k].node_count_max, a, b, h, eps );
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

