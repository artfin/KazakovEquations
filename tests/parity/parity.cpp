#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

enum class Parity
{
    EMPTY,
    ODD,
    EVEN
};

std::ostream & operator<<( std::ostream & os, const Parity & p )
{
    if ( p == Parity::EVEN )
        os << "EVEN";
    else if ( p == Parity::ODD )
        os << "ODD";
    else if ( p == Parity::EMPTY )
        os << "NO PARITY";

    return os;
}

std::istream & operator>>( std::istream & os, Parity & p )
{
    std::string str;
    std::stringstream ss;
    
    ss << os.rdbuf();
    str = ss.str();

    str.erase(str.find_last_not_of(" \n\r\t")+1);    
    std::cout << "(operator>>) str = " << str << std::endl;

    if ( str == "EVEN" )
    {
        std::cout << "setting parity to EVEN" << std::endl;
        p = Parity::EVEN;
    }
    else if ( str == "ODD" )
        p = Parity::ODD;
    else if ( str == "NO PARITY" )
        p = Parity::EMPTY;

    return os;
}

int main()
{
    std::ifstream ifs( "./db.txt" );

    Parity p;
    ifs >> p;

    std::cout << "parity: " << p << std::endl;

    return 0;
}

