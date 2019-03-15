#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>

class Eigenvalue
{
public:
    Eigenvalue( )
    {
    }

    Eigenvalue( double value, int J, int M, int node_count ) :
        value(value), J(J), M(M), node_count(node_count)
    {
    }

    bool operator<( const Eigenvalue & other ) const
    {
        return (this->J < other.getJ()) || 
              ((this->J == other.getJ()) && (this->M < other.getM())) ||
              ((this->J == other.getJ()) && (this->M == other.getM()) && (this->node_count < other.getNodeCount()));        
    }

    friend std::ifstream& operator>>( std::ifstream & os, Eigenvalue & eig ); 
    friend std::ostream& operator<<( std::ostream & os, const Eigenvalue & eig );

    void setAngularMomentum( int J, int M )
    {
        this->J = J;
        this->M = M;
    }

    void setValue( double value ) { this->value = value; }
    void setNodeCount( int node_count ) { this->node_count = node_count; }

    double getValue() const { return value; }
    int getM() const { return M; }
    int getJ() const { return J; }
    int getNodeCount() const { return node_count; }

private:
    double value;
    int J;
    int M;

    int node_count;
};

std::ifstream& operator>>( std::ifstream & os, Eigenvalue & eig )
{
    double value;
    int J, M, node_count;
    os >> J >> M >> node_count >> value;

    eig.setAngularMomentum( J, M );
    eig.setValue( value );
    eig.setNodeCount( node_count );

    return os;
} 

std::ostream& operator<<( std::ostream & os, const Eigenvalue & eig )
{
    os << eig.getJ() << " " << eig.getM() << " " << eig.getNodeCount() << " " << eig.getValue(); 
    return os;
}

class EigenvalueDatabase
{
public:
    EigenvalueDatabase( std::string db_filename ) : db_filename(db_filename)
    {
    }

    void read( )
    {
        if ( is_file_exists(db_filename) )
        {
            std::ifstream ifs( db_filename );

            for ( Eigenvalue tmp; ifs >> tmp; db.push_back(tmp) );
        }
        else
        {
            std::cerr << "Database doesn't yet exist! Nothing to read." << std::endl;
        } 
    }

    bool is_file_exists( std::string filename )
    {
        std::ifstream infile( filename );
        return infile.good();
    }

    bool check_unique( const Eigenvalue & eig )
    // checks that there is no value in database with 
    // the same exact {J, M, node_count}
    {
        for ( Eigenvalue & eig_db : db )
        {
            if ( (eig_db.getJ() == eig.getJ()) &&
                 (eig_db.getM() == eig.getM()) && 
                 (eig_db.getNodeCount() == eig.getNodeCount()) )
                return false; 
        }

        return true;
    }

    void show()
    {
        std::cout << "-----------------------" << std::endl;
        for ( Eigenvalue const & eig : db )
           std::cout << eig << std::endl; 

        std::cout << "Size: " << db.size() << std::endl;
        std::cout << "-----------------------" << std::endl;
    }

    void insert( std::vector<Eigenvalue> const & eigs )
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

    void write()
    {
        std::cout << "Writing eigenvalues to the DB." << std::endl;
        
        std::sort( db.begin(), db.end() );

        std::ofstream ofs( db_filename );

        for ( Eigenvalue & eig : db )
            ofs << eig << std::endl;
    }


private:
    std::string db_filename; 

    std::vector<Eigenvalue> db;
};

int main()
{
    std::string db_filename = "./eigenvalues.db";

    EigenvalueDatabase EigenDB( db_filename );
    EigenDB.read();

    EigenDB.show();

    Eigenvalue e1( -0.00554, 6, 0, 3 );
    Eigenvalue e2( -0.00375, 5, 0, 2 ); 
    Eigenvalue e3( -0.00375, 5, 0, 1 ); 
    Eigenvalue e4( -0.00375, 5, 1, 1 ); 
    Eigenvalue e5( -0.00375, 5, 1, 2 ); 
    Eigenvalue e6( -0.00375, 4, -1, 3 ); 
    Eigenvalue e7( -0.00375, 4, 0, 3 ); 
    std::vector<Eigenvalue> eigs{ e1, e2, e3, e4, e5, e6, e7 };

    EigenDB.insert( eigs ); 

    EigenDB.write();

    return 0;
}
