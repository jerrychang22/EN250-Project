#ifndef _STRUCTS_HPP_
#define _STRUCTS_HPP_

#include "simulation.hpp"
#include "species.hpp"

//#include "constants.hpp"
#include <memory>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <exception>
#include <memory>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iterator>


#include "output.hpp"

//class World {
//  // this class contains all the variables/switches of the simulation environment:
//  // - list of species
//  // - list of trees
//  // - switch: how is the ligth treated (from above / from somewhere else?)
//  // - switch: is the ground flat?
//  // also, this class contains initialization and loading/saving procedures
//  public:
//
//};




///*======*///
///* TREE *///
///*======*///
class Tree {
  public:
    // useful for ulterior reference purposes:
    Simulation* S;
    /* tree index */
    unsigned int n;
    /* species index */
    double s;
    /* species pointer */
    //Species* species;
    std::shared_ptr <Species> species;
    /* age */
    int age;
    /* understory condition: 0 := canopy, 1 := understory, 2 := area in both */
    unsigned short is_ustory;
    /* root position */
    unsigned short xo;
    unsigned short yo;
    /* center of growth position */
    unsigned short x;
    unsigned short y;
    //    unsigned short z; /* altitude (jean: not used - will be ever useful ?) */
    /* short dx[Sim.SECTORS];*/ /* farthest distance per sector */
    /* short dy[Sim.SECTORS];*/ /* farthest distance per sector */
    double d; /* trunk diameter (dm) */
    double h; /* tree height (dm) */
    double w; /* crown radius (dm) */
    double cbh; /* crown base height */

    int SECTORS; // the number of sectors of that tree

    std::vector<double> ws; /* crown sector radii (dm) */
    double max_ws; /* (Temporary Variable) max crown radius across sectors (dm) */
    std::vector<double> ustory_ws; /* understory sector radii (dm) */
    std::vector<double> max_sectors; /* used to store max value of over/understory */
    std::vector<double> cbhs; /* crown sector base heights */
    std::vector<unsigned short> sector_story; /* indicator of whether this sector is 0,1,2 */
    /* double gr; */ /* growth rate  */ /* constant growth rate for flat top model */
    double mr; /* mortality rate */
    double tr; /* crown transparency */
    double pot_ca; /* potential crown area calculated in pre_crown */
    unsigned int ca; /* crown area (dm2 = No of points) */
    unsigned int ustory_ca; /* understory crown area (dm2 = No of points) */
    unsigned int brdl; /* length of the crown border */
    double brdh; /* average height of the crown border */

    double pot_radius;

    double root_w; // the root radius (dm)
    double water_uptake=1.0; // the water uptake
    double water_needed=1.0; // the water uptake needed

#ifdef __USE_ROOTS
    double wdemand; /* water demand of tree */
    double wuptake; /* the most water a tree can use */
    /*v the realized water uptake; 0: top, 1: mid, 2: bot, 3: total for tree */
    double wrealuptake[ROOTMAX]; /* Convenient that 3 layers of roots */
    long double watersect[Sim.SECTORS]; /* realized sector uptake top level */
    long double watersectmid[Sim.SECTORS]; /* realized sector uptake mid level */
    long double watersectbot[Sim.SECTORS]; /* realized sector uptake bot level */
    short wlimit; /* 0 := light limited, 1 := water limited */
#endif
    struct PIXEL color;  /* canopy color on biomass graph */

    double** claims;

    double ComputeClaim(double r2, double w2, int sector)
    {
      // dist is a squared ratio: [claim radius/crown radius]**2
      if ( S->RESTRICTCROWN == 0 ) {
        // the shape starts from the ground
        return (get_claim_new(r2/w2, species->csexp, h, species->gap));
      } else if ( S->RESTRICTCROWN == 1 ) {
        // the shape starts from tree[t]->cbhs[sector]
        return (cbhs[sector] + get_claim_new(r2/w2, species->csexp, h-cbhs[sector], species->gap));
      } else if ( S->RESTRICTCROWN == 2 ) {
        // the shape starts from the ground, and is then cropped at tree[t]->cbhs[sector]
        double v = get_claim_new(r2/w2, species->csexp, h, species->gap);
        if (v < cbhs[sector]) {
          return (-999); // arbitrary high negative height to cancel the claim
        }
      } else if ( S->RESTRICTCROWN == 3 ) {
        // the shape starts from the ground, and is computed according to the max possible growth
        //double d_ratio = d/species->dMax;
        //if (d_ratio > 1) { d_ratio = 1; }
        //return (get_claim_new(r2/ (exp(log(species->wMax * d_ratio)*2)), species->csexp, h, species->gap));
        //return (get_claim_new(r2/ (exp(log(species->wMax)*2)), species->csexp, h, species->gap));
        double max_max_sectors = 0;
        for (int i=0;i<SECTORS;i++){ if (max_sectors[i] > max_max_sectors) {max_max_sectors = max_sectors[i];}}
        return (get_claim_new(r2/ (max_max_sectors * max_max_sectors), species->csexp, h, species->gap));
      }
      return -1;
    }

    double get_claim_new(double r2_w2, double e, double h, double gap) { /* compute the claim */
      double shape;
      if (S->CROWNSHAPE == 0) {
        /* CASE 0: h = h_max * (1 - (r / r_max)^exp) */
        shape = 1.0 - exp(log(r2_w2) * e / 2.0);
      } else {
        /* CASE 1: h = h_max * (1 - r / r_max)^(1 / exp) */
        if (e>=100) {
          shape = 1.0;
        } else {
          shape = exp(log(1.0 - sqrt(r2_w2)) / e);
        }
      }
      return h * (gap + (1.0 - gap) * shape);
    }



    //// Variant with unique ownership
    //static void CreateTree(World& world)
    //{
    //  std::unique_ptr<Tree> p(new Tree());
    //  //world.trees.push_back(std::move(p));
    //}

    Tree() {}
    Tree(Simulation* Si) : ws(Si->SECTORS), ustory_ws(Si->SECTORS), max_sectors(Si->SECTORS), cbhs(Si->SECTORS), sector_story(Si->SECTORS) {S = Si; SECTORS = Si->SECTORS;}

    //private:


    //private:

    //  // Declare *structors private to avoid possibility to create Entity instances
    //  // without CreateGameEntity() method, e.g. on stack.
    //  Tree();
    //  Tree(const Tree&);
};

#ifdef __USE_ROOTS
/* ============== */
/* ROOT PARAMTERS */
/* ============== */
struct GROUND { /* root data */
  /* root index (same as tree index) */
  unsigned int n;
  /* species index (same as tree index) */
  unsigned short s;
  /* main position (will not change, once assigned) */
  unsigned short x;
  unsigned short y;
  /* sector root areas (dm2 = No of points) */
  unsigned int rarea_top[Sim.SECTORS];
  unsigned int rarea_mid[Sim.SECTORS];
  unsigned int rarea_bot[Sim.SECTORS];
  /* total areas (sum of sector areas) */
  unsigned int total_top;
  unsigned int total_mid;
  unsigned int total_bot;
  /* top root sector radii (dm) */
  double radii_top[Sim.SECTORS];
  double radii_mid[Sim.SECTORS];
  double radii_bot[Sim.SECTORS];
  /* maximum radius for "super root area" (dm) */
  double radius;
};
#endif /* __USE_ROOTS */

#endif
