#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include <cmath> // to define M_PI
#include <string> // to access std::string
#include "species.hpp" // we have a vector of species

//#include <vector> // std::vector
#include <deque> // std::deque
#include <memory> // std::unique_ptr


class Simulation {
  public:


    /* ========================== */
    /* GENERAL SIMULATION OPTIONS */
    /* ========================== */

    int HSCALE; /* height profile graph scale (dm) NO LESS THAN THE HIGHEST SPECIES */
    int CYCLES; /* simulation length */
    int SEED; /* random seed - use 0 to switch to time(NULL) */
    int SECTORS; /* number of sectors per tree: 0: off; 4, 8: on */
    int SMAX; /* number of species in Model */
    double SEEDS_PER_HA; // yearly seedling density


    /* =================== */
    /* ENVIRONMENT OPTIONS */
    /* =================== */

    int XMAX; /* should be 2^n */
    int YMAX; /* should be 2^n */
    int ZMAX; /* (in dm) should be at least HSCALE, perhaps more if there is a slope */
    double WATER_AVAILABLE; /* (in dm) should be at least HSCALE, perhaps more if there is a slope */

    // COMPLETEDISTURBANCE
    // value from 0 to 1, corresponds to the annual probability of complete disturbance
    double COMPLETEDISTURBANCE;

    /* =================== */
    /* SOLAR MODEL OPTIONS */
    /* =================== */

    double CORIENTATION; /* Slope orientation: highest point angle, in radians. (0.5*M_PI): North, (M_PI): West, (1.5*M_PI): South, (0): East */
    double CANGLE; /* Slope's angle (in radians, should be any value in [0,0.5*M_PI[) - 0 indicates a flat ground */
    double CSUNALT; /* Sun altitude (= elevation) in radians, should be any value in ]0,0.5*M_PI] - 0 indicates a sun hidden by the horizon, 0.5*M_PI indicates zenith */
    double CSUNORIENTATION; /* Sun azimuth. Note that this value is NOT as in the standard horizontal coordinate system, here we measure it in radians and counter clockwise. (0.5*M_PI): (light coming from) North, (M_PI): West, (1.5*M_PI): South, (0): East */
    // for crater lake (lat=42.9388312° N, long=122.1459961° W):
    // the slope angle is approx 36 degree (max altitude of the cone is 755 feet, with approx 1000 feet radius),   so ANGLE (0.2*M_PI)
    // on 21/09/2013, at 13:00 PM, elevation=47.46° (slightly more than Pi/4), azimuth=179.45° (South)
    // => zenith is for SUNORIENTATION (1.54*M_PI)   and   SUNALT (0.26*M_PI)


    /* ================== */
    /* MECHANISMS OPTIONS */
    /* ================== */

    // ROOT
    // 0 := Off, no root system at all
    // 1 := explicit root system
    // 2 := implicit root system
    int ROOT;
    //int ROOTMAX 4 /* maximum number of trees that can have roots claim a spot */
    
    // LIGHTWATERGROWTH
    // 0 := only light limitation
    // 1 := multiplicative influence (taking the product (light limitation) * (water limitation))
    // 2 := Liebig's law of the minimum
    int LIGHTWATERGROWTH;

    // LIGHTWATERMORTALITY
    // 0 := only light limitation (position relative to canopy)
    // 1 := additive influence (taking the sum of light limitation and water limitation)
    // 2 := Liebig's law of the minimum (here, "maximum")
    int LIGHTWATERMORTALITY;

    // RESTRICTCROWN
    // 0 := Off, the crown base height (tree[t].cbhs[i]) has no influence on the shape
    // 1 := the shape starts at tree[t].cbhs[i]. This will lead to a tree with a completely flat arborescence if the competition is enough to force the crown base height to "climb" to the top.
    // 2 := the shape is computed from the ground, but is cropped at tree[t].cbhs[i].
    // 3 := the shape is computed from the ground, according to the optimal diameter/crown ratio.
    int RESTRICTCROWN; 
    // !! RESTRICTCROWN == 2 is not performing well anymore (cbhs is initialized at 2.00 in the beginning of grow_up())
    // !! RESTRICTCROWN == 3 is not performing well (if all sectors are computed according to the biggest, then they are mostly flat on top -> leading to weird shapes)

    // VERTICALGROWTH
    // 0 := vertical growth is dependent on trunk diameter
    // 1 := vertical growth is dependent on the vertical proportion of the sector that is sunlit
    int VERTICALGROWTH;

    // RESTRICTHEIGHT
    // this constant is used only when VERTICALGROWTH != 0
    // 0 := height increment is constant (no restriction, so tree[t].h[s] can be more than species[i].hMax)
    // 1 := height increment is constant until the tree height reaches species[i].hMax (i.e. capped)
    // 2 := height increment is proportional to | species[i].hMax - current height |
    int RESTRICTHEIGHT;

    // RESTRICTRADIUS
    // 0 := radius increment is constant (no restriction, so tree[t].ws[s] can be more than species[i].wMax)
    // 1 := radius increment is constant until the sector radius reaches species[i].wMax (i.e. capped)
    // 2 := radius increment is proportional to | species[i].wMax - current radius at crown base |
    // 3 := radius increment is proportional to | maximal width allowed at crown base height - current radius at crown base |
    //     note that 3 is equivalent to 2 when RESTRICTCROWN == 1
    int RESTRICTRADIUS;

    // POTENTIALCROWN
    // 0 := the potential crown is computed based on the crown size from each sectors at last timestep (this leads to very small potential crowns)
    // 1 := the potential crown is computed based on the crown size from the max sector at last timestep (logic behind: tree structure must be suitable to accomodate its largest crown extent)
    // 2 := the potential crown is increased at each timestep, by considering that the tree is not in competition (idealized potential crown)
    int POTENTIALCROWN;

    int SEEDANYWHERE; /* 0 := Use newLight & newShade parameters, 1 := Seed anywhere */
    int INTERMEDIATESPECIES; /* 0 := default, use pre-defined species, 1 := generate a random intermediate species for each tree */
    int CROWNSHAPE; /* 0 := h = h_max * (1 - r^exp / r_max^exp)[Drew's profile if exp == 2]; 1 := h = h_max * (1 - r / r_max)^(1/exp)[SORTIE profile if exp == 2] */
    //int VORONOI = 0; /* 0 := Off, 1 := On */
    int REAL3D; /* 0 := Off, 1 := On */


    /* ============== */
    /* BITMAP OPTIONS */
    /* ============== */

    int OUTPUTMOD; /* 0: output all iterations; n>1: output every n iterations */
    int OUTPUTBMP; /* 0 := False, 1 := True */
    int OUTPUTALLTREES; /* 0 := False, 1 := True */
    int OUTPUTSTATS; /* 0 := False, 1 := True */
    int SHOWTREELIMITS; /* 0 := Hide, 1 := Show */
    int SHOWTREEROOT; /* 0 := Hide, 1 := Show */
    int SHOWTREECENTER; /* 0 := Hide, 1 := Show */


    /* ========================= */
    /* PROGRAM EXECUTION OPTIONS */
    /* ========================= */

    int REDO_XP; /* 0 := ask, 1 := redo, 2 := do nothing, 3 := redo only if Stats.csv cannot be found or contains less than CYCLES lines */
    std::string output_dirname;

    /* =============================== */
    /* OTHERS (initialized at runtime) */
    /* =============================== */

    double nSeedTotal;
    //std::vector <Species> species; // the vector containing the different species
    //std::deque <Species> species; // the vector containing the different species
    std::deque<std::shared_ptr<Species>> species;


    /* ========= */
    /* GRAVEYARD */
    /* ========= */

    //#ifdef __USE_ROOTS
    //      int SHOWWATER 1 /* 0 := Hide, 1:= Show */
    //#else
    //      int SHOWWATER 0
    //#endif

    // ZRnd
    // double value allowing a tree to claim with 50% chance a tree which is ZRnd dm higher than himself
    // set at 0.0 to deactivate
    // double ZRnd;
    // NOT USED ANYMORE


    /* ======= */
    /* METHODS */
    /* ======= */

    double initialize_nSeedTotal();


    /* ============ */
    /* CONSTRUCTORS */
    /* ============ */
    
    
    //Simulation();
    Simulation(int ac, char** av);

};

#endif
