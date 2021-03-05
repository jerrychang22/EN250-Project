#ifndef _SPECIES_HPP_
#define _SPECIES_HPP_

#include "output.hpp" // to use PIXEL
#include <vector> // std::vector

/* ======= */
/* SPECIES */
/* ======= */
class Species {
  // this class gives all parameters of a species
  public:
    double hMax; /* maximum height (dm) */
    double dMax; /* maximum trunk diameter (dm) */
    double wMax; /* maximum crown radius (dm) */
    double gap;  /* fraction of height free of crown */
    double dInit; /* initial trunk diameter (dm) */
    double tiltMax; /* tangent of the maximum trunk canting angle */
    /* double gRate; */ /* growth rate */
    double mrLight; /* mortality rate for trees with canopy */
    double mrShade; /* mortality rate for trees without canopy */
    double newLight; /* probability to germinate in open space */
    double newShade; /* probability to germinate in shade */
    double sunShade; /* fraction of light in the shade */
    /* PRI parameters */
    double b1;
    double b2;
    double b3;
    double csexp; /* crown shape exponent (1 for conifer, 2 for round) */
    double rinc; /* yearly radius increment (species specific) [crown x-y area] */
    double hinc; /* yearly height increment (species specific) [crown z distance] */  // plausible order of magnitude, see figs. 100 and 101, p. 199 of TPOFYS
    double hd; // TODO this is an adhoc way of boosting height in shade intolerant species
    double water_mr;
    double water_needed;
    double shadetolerance;
    double growth_factor;
    //double rootinc; /* yearly root increment (species specific) */
    //double waterdep[3]; /* how dependent upon water the species is */
    double nSeed; /* number of new trees planted */
    struct PIXEL leaf; /* height map leaf color */
//#if SHOWWATER
//    struct PIXEL waterleaf;
//#endif
    // this class gives seed dispersal procedures too
    // private:

    Species(int s); // generate a pre-defined species (default to Eastern Hemlock when s==0, and White Pine when s==1)
    Species(double alpha); // generate an intermediate species with parameters between Eastern Hemlock (s==0) and White Pine (s==1)
};

#endif
