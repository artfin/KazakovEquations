#pragma once

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

std::ostream& operator<<( std::ostream& os, const Parity & p );
std::istream& operator>>( std::istream& os, Parity & p );

