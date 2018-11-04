#include "./parity.hpp"

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

std::istream & operator>>( std::istream & is, Parity & p )
{
    std::string str;
    is >> str;

    if ( str == "EVEN" )
        p = Parity::EVEN;
    else if ( str == "ODD" )
        p = Parity::ODD;
    else if ( str == "NO PARITY" )
        p = Parity::EMPTY;

    return is;
}
