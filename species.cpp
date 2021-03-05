#include "species.hpp"
#include "simulation.hpp"


#include <iostream> // for debug msg

//#include <cmath>

Species::Species(int s)
{
  hMax = 250; /* maximum height (dm) */
  dMax = 10; /* maximum trunk diameter (dm) */
  //  wMax = 75; /* maximum crown radius (dm) */
  gap = 0.1;  /* fraction of height free of crown */
  dInit = 0.4; /* initial trunk diameter (dm) */
  tiltMax = 0.01; /* tangent of the maximum trunk canting angle */
  mrLight = 0.01; /* mortality rate for trees with canopy; */
  newLight = 1; /* probability to germinate in open space */
  newShade = 0.01; /* probability to germinate in shade */
  sunShade = 0.09; /* fraction of light in the shade */
  csexp = 2; /* crown shape exponent (1 for conifer, 2 for round) */
  // rinc = 2; /* yearly radius increment (species specific) [crown x-y area] */
  rinc = 15; /* yearly radius increment crown factor */
  hinc = 5; /* yearly height increment (species specific) [crown z distance] */  // plausible order of magnitude, see figs. 100 and 101, p. 199 of TPOFYS
  water_mr = 0; // mortality factor induced by water limitation
  water_needed = 1; // need of a tree per surface unit of roots
  shadetolerance = 1; // by default the species is shade tolerant
  growth_factor = 1.; // no alteration of growth by default
  if (s==0) {
    // Species 0 is by default the shade tolerant species: slow growth, low mortality (under shade), low max height
    // parameterization is suited for Eastern Hemlock
    leaf.R = 0xB0; // green
    leaf.G = 0xFF;
    leaf.B = 0x00;
    b1 = 1.749000; // b1 = b1 * exp(log(10) * b2) / 5;
    b2 = -0.436216; // b2 = b2 + 1;
    b3 = 0.968043; // b3 = exp(log(b3) * 10);
    mrShade = 0.01; /* mortality rate for trees without canopy; */
    shadetolerance = 1;
    wMax = 126; // 12m60 max crown radius, after 150 years of full growth
    dMax = 8.7734619; // 0.88 meters, after 150 years of full growth
  } else if (s==1) {
    // Species 1 is by default the shade intolerant species: fast growth, high mortality (under shade), high max height
    // parameterization is suited for White Pine
    leaf.R = 0xFF; // yellow
    leaf.G = 0xB0;
    leaf.B = 0x00;
    b1 = 3.856766; // b1 = b1 * exp(log(10) * b2) / 5;
    b2 = -0.553002; // b2 = b2 + 1;
    b3 = 0.979675; // b3 = exp(log(b3) * 10);
    mrShade = 0.2; /* mortality rate for trees without canopy; */
    shadetolerance = 0;
    wMax = 208; // 20m80 max crown radius, after 150 years of full growth
    dMax = 14.2793204; // 1.43 meters, after 150 years of full growth
  } else {
    // we leave to the configuration file the task to initialize variables
    leaf.R = rand() % 0x100;
    leaf.G = rand() % 0x100;
    leaf.B = rand() % 0x100;
  }
  b1 = b1 * exp(log(10.) * b2) / 5.;
  b2 = b2 + 1.;
  b3 = exp(log(b3) * 10.);
}


double generate_intermediate()
{
  // inspired by the polynomial mutation operator (Deb and Goyal 1996)
  double rnd = (double) rand() / RAND_MAX;
  double y;
  if ((double) rand() / RAND_MAX < 0.5) {
    y = 0.05;
  } else {
    y = 0.95;
  }
  if ((double) rand() / RAND_MAX < 0.5) {
    y = 0.;
  } else {
    y = 1.;
  }
  double yl = 0.;
  double yu = 1.;
  double delta1 = (y-yl)/(yu-yl);
  double delta2 = (yu-y)/(yu-yl);
  double eta_m = 20;
  double mut_pow = 1.0/(eta_m+1.0);
  double xy,val,deltaq;
  if (rnd <= 0.5)
  {
    xy = 1.0-delta1;
    val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
    deltaq =  pow(val,mut_pow) - 1.0;
  }
  else
  {
    xy = 1.0-delta2;
    val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
    deltaq = 1.0 - (pow(val,mut_pow));
  }
  y = y + deltaq*(yu-yl);
  if (y<yl)
    y = yl;
  if (y>yu)
    y = yu;
  return(y);
}


Species::Species(double alpha)
{
  // this generates a new species, somewhere between prototypical species 1 (Eastern Hemlock) and species 2 (white pine)
  // all parameters are copy-n-paste from the other function def
  double alpha1, alpha2;
  if (alpha <= 1.) {
    alpha1 = alpha;
    alpha2 = alpha;
  } else if (alpha <= 2.) {
    alpha1 = alpha - 1.;
    if (alpha < 1.5) {
      alpha2 = 0.;
    } else {
      alpha2 = 1.;
    }
  } else {
    alpha1 = generate_intermediate();
    alpha2 = alpha1;
    // std::cout << "just created an intermediate species with polynomial mutation and y=" << alpha1 << std::endl;
  }
  hMax = 250;
  dMax = 10;
  gap = 0.1;
  dInit = 0.4;
  tiltMax = 0.01;
  mrLight = 0.01;
  newLight = 1;
  newShade = 0.01;
  sunShade = 0.09;
  csexp = 2;
  rinc = 15;
  hinc = 5;
  water_mr = 0;
  water_needed = 1;
  shadetolerance=1;

  leaf.R = 0xB0 + (unsigned char)(alpha1*(0xFF - 0xB0)); // s==1 => this is White Pine; s==0 => this is Eastern Hemlock
  leaf.G = 0xFF - (unsigned char)(alpha1*(0xFF - 0XB0));
  leaf.B = 0x00;
  b1 = 1.749000 + alpha1*(3.856766 - 1.749000);
  b2 = -0.436216 + alpha2*(0.436216 - 0.553002);
  b3 = 0.968043 + alpha2*(0.979675 - 0.968043);
  mrShade = 0.01 + alpha1*(0.2 - 0.01);
  shadetolerance = 1 - alpha1;
  wMax = 126. + alpha2*(208. - 126.);
  dMax = 8.7734619 + alpha2*(14.2793204 - 8.7734619);

  b1 = b1 * exp(log(10) * b2) / 5;
  b2 = b2 + 1;
  b3 = exp(log(b3) * 10);
}



//Species::Species(Simulation* Si)
//{
//    leaf.R = 0xB0;
//    leaf.G = 0xFF;
//    leaf.B = 0x00;
//    nSeed = Sim->nSeedTotal / Sim->SMAX;
//    b1 = b1 * exp(log(10) * b2) / 5;
//    b2 = b2 + 1;
//    b3 = exp(log(b3) * 10);
//}

