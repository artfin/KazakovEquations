#pragma once

namespace constants
{
    // convert Hartree to inverse centimeters 
    //const double HTOCM = 219474.6313702; 
    
    // bohr to angstrom
    //const double BOHRTOANG = 0.52917721067;

    // Hartree to inverse cm
    const double HTOCM = 219475.797;

    // atomic length unit to meters == a0
    const double ALU = 0.52917721067 * 1E-10;

    // bohr to angstrom
    const double BOHRTOANG = 0.529177;
    const double alpha_inv = 137.035999139;
}

namespace dunker 
{
    const double P0R = 1.000;
    const double P1R = 0.510;
    const double P2R = 0.780;

    const double P0A = 1.000;
    const double P1A = 0.320;
    const double P2A = 0.240;

    //double epsilon = 140.39; // cm-1 ?
    const double epsilon = 140.39566; // cm-1 ?
    const double Rm = 3.930; // angstroms
    const double alpha = 13.5;

}
