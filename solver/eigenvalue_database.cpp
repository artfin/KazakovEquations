#include "./eigenvalue_database.hpp"

bool EigenvalueDatabase::file_exists( std::string const & filename )
{
    std::ifstream infile( filename );
    return infile.good();
}

void EigenvalueDatabase::read( )
{
    if ( file_exists(db_filename) )
    {
        std::ifstream ifs( db_filename );
        for ( Eigenvalue tmp; ifs >> tmp; db.push_back(tmp) );
        //{
            //std::cout << tmp << std::endl;
        //}
    }
    else
    {
        std::cerr << "Database doesn't yet exist! Nothing to read." << std::endl;
    } 
}

bool EigenvalueDatabase::check_unique( const Eigenvalue & eig )
// checks that there is no value in database with 
// the same exact {J, M, node_count}
{
    for ( Eigenvalue & eig_db : db )
    {
        if ( (eig_db.get_J() == eig.get_J()) &&
                (eig_db.get_M() == eig.get_M()) && 
                (eig_db.get_node_count_min() == eig.get_node_count_min()) )
            return false; 
    }

    return true;
}

void EigenvalueDatabase::show()
{
    std::cout << "-----------------------" << std::endl;
    for ( Eigenvalue const & eig : db )
        std::cout << eig << std::endl; 

    std::cout << "Size: " << db.size() << std::endl;
    std::cout << "-----------------------" << std::endl;
}

void EigenvalueDatabase::insert( std::vector<Eigenvalue> const & eigs )
{
    for ( Eigenvalue const & eig : eigs )
    {
        bool unique = check_unique( eig );   
        if ( !unique )
        {
            std::cerr << "There is already an eigenvalue with the same set of numbers in database!" << std::endl;    
            std::cerr << "Eigenvalue: " << eig << std::endl;
        } 
        else
        {
            db.push_back( eig ); 
        }
    }
}

void EigenvalueDatabase::write()
{
    std::cout << "Writing eigenvalues to the DB." << std::endl;

    std::sort( db.begin(), db.end() );

    std::ofstream ofs( db_filename );
    ofs << std::fixed << std::setprecision(13);
        
    for ( Eigenvalue & eig : db )
        ofs << eig << std::endl;
}

