#ifndef _CONSTANTS_HPP_
#define _CONSTANTS_HPP_

#define ANGLESIM 15.0





// the following switches are used for output of the ESA meeting poster
#define ROOT
//#define WATER_GRADIENT
//#define HIGHER_DIST
//#define LOWER_DIST

//#define ESA_WATER
//#define LOW_WATER
//#define MID_WATER
//#define HIGH_WATER
//#define OTHER_WATER





/* ================ */
/* OUTPUT CONSTANTS */
/* ================ */

// #define DEBUG
#ifdef DEBUG
 #define D 
#else
 #define D for(;0;)
#endif


/* ==================== */
/* SIMULATION CONSTANTS */
/* ==================== */

// // // #define HSCALE 250 /* height profile graph scale (dm) NO LESS THAN THE HIGHEST SPECIES */
// // // #define CYCLES 1500 /* simulation length */

//#define XMAX 512 /* (in dm) should be 2^n */
//#define YMAX 2048 /* (in dm) should be 2^n */
//#define XMAX 1024 /* (in dm) should be 2^n */
//#define YMAX 1024 /* (in dm) should be 2^n */
// // // #define XMAX 512 /* should be 2^n */
// // // #define YMAX 512 /* should be 2^n */
//#define XMAX 64 /* should be 2^n */
//#define YMAX 64 /* should be 2^n */

// // // #define ZMAX 512 /* (in dm) should be at least HSCALE, perhaps more if there is a slope */
// // // #define TMAX 7000 /* maximum number of trees */
// // // #define SECTORS 8 /* number of sectors per tree: 0: off; 4, 8: on */
//#define CORIENTATION ( ANGLESIM *2.0*M_PI /16.0) /* Slope orientation: highest point angle, in radians. (0.5*M_PI): North, (M_PI): West, (1.5*M_PI): South, (0): East */
// // // #define CORIENTATION (0.5*M_PI) /* Slope orientation: highest point angle, in radians. (0.5*M_PI): North, (M_PI): West, (1.5*M_PI): South, (0): East */

//#define CANGLE (0.2*M_PI) /* Slope's angle (in radians, should be any value in [0,0.5*M_PI[) - 0 indicates a flat ground */
// // // #define CANGLE (0.) /* Slope's angle (in radians, should be any value in [0,0.5*M_PI[) - 0 indicates a flat ground */

//#define CSUNALT (0.26*M_PI) /* Sun altitude (= elevation) in radians, should be any value in ]0,0.5*M_PI] - 0 indicates a sun hidden by the horizon, 0.5*M_PI indicates zenith */
// // // #define CSUNALT (0.5*M_PI) /* Sun altitude (= elevation) in radians, should be any value in ]0,0.5*M_PI] - 0 indicates a sun hidden by the horizon, 0.5*M_PI indicates zenith */

// // // #define CSUNORIENTATION (1.54*M_PI) /* Sun azimuth. Note that this value is NOT as in the standard horizontal coordinate system, here we measure it in radians and counter clockwise. (0.5*M_PI): (light coming from) North, (M_PI): West, (1.5*M_PI): South, (0): East */
//#define SUNORIENTATION (1.54*M_PI) /* Sun azimuth. Note that this value is NOT as in the standard horizontal coordinate system, here we measure it in radians and counter clockwise. (0.5*M_PI): (light coming from) North, (M_PI): West, (1.5*M_PI): South, (0): East */
// for crater lake (lat=42.9388312° N, long=122.1459961° W):
// the slope angle is approx 36 degree (max altitude of the cone is 755 feet, with approx 1000 feet radius),   so ANGLE (0.2*M_PI)
// on 21/09/2013, at 13:00 PM, elevation=47.46° (slightly more than Pi/4), azimuth=179.45° (South)
// => zenith is for SUNORIENTATION (1.54*M_PI)   and   SUNALT (0.26*M_PI)


/* ============= */
/* MODEL OPTIONS */
/* ============= */

// // // #define SMAX 2 /* number of species in Model */

// // // #define ROOTMAX 4 /* maximum number of trees that can have roots claim a spot */

// RESTRICTCROWN
// 0 := Off, the crown base height (tree[t].cbhs[i]) has no influence on the shape
// 1 := the shape starts at tree[t].cbhs[i]. This will lead to a tree with a completely flat arborescence if the competition is enough to force the crown base height to "climb" to the top.
// 2 := the shape is computed from the ground, but is cropped at tree[t].cbhs[i].
// // // #define RESTRICTCROWN 0 
// !! RESTRICTCROWN == 2 is not performing well anymore (cbhs is initialized at 2.00 in the beginning of grow_up())

// VERTICALGROWTH
// 0 := vertical growth is dependent on trunk diameter
// 1 := vertical growth is dependent on the vertical proportion of the sector that is sunlit
// // // #define VERTICALGROWTH 0

// RESTRICTHEIGHT
// this constant is used only when VERTICALGROWTH != 0
// 0 := height increment is constant (no restriction, so tree[t].h[s] can be more than species[i].hMax)
// 1 := height increment is constant until the tree height reaches species[i].hMax (i.e. capped)
// 2 := height increment is proportional to | species[i].hMax - current height |
// // // #define RESTRICTHEIGHT 2

// RESTRICTRADIUS
// 0 := radius increment is constant (no restriction, so tree[t].ws[s] can be more than species[i].wMax)
// 1 := radius increment is constant until the sector radius reaches species[i].wMax (i.e. capped)
// 2 := radius increment is proportional to | species[i].wMax - current radius at crown base |
// 3 := radius increment is proportional to | maximal width allowed at crown base height - current radius at crown base |
//     note that 3 is equivalent to 2 when RESTRICTCROWN == 1
// // // #define RESTRICTRADIUS 1

// ZRnd
// double value allowing a tree to claim with 50% chance a tree which is ZRnd dm higher than himself
// set at 0.0 to deactivate
// // // #define ZRnd 1.0


// // // #define SEEDANYWHERE 1 /* 0 := Use newLight & newShade parameters, 1 := Seed anywhere */
// // // #define CROWNSHAPE 1 /* 0 := h = h_max * (1 - r^exp / r_max^exp)[Drew's profile if exp == 2]; 1 := h = h_max * (1 - r / r_max)^(1/exp)[SORTIE profile if exp == 2] */
// // // #define VORONOI 0 /* 0 := Off, 1 := On */
// // // #define REAL3D 0 /* 0 := Off, 1 := On */


/* ================ */
/* BITMAP CONSTANTS */
/* ================ */

// // // #define OUTPUTBMP 1 /* 0 := False, 1 := True */
// // // #define SHOWTREELIMITS 0 /* 0 := Hide, 1 := Show */
// // // #define SHOWTREEROOT 1 /* 0 := Hide, 1 := Show */
// // // #define SHOWTREECENTER 1 /* 0 := Hide, 1 := Show */
#ifdef __USE_ROOTS
#define SHOWWATER 1 /* 0 := Hide, 1:= Show */
#else
#define SHOWWATER 0
#endif


//#define __USE_ROOTS






// REWRITE CPP => this function is still used in mod.cpp (eventually it shouldnt) and in structs.hpp (eventually it should be moved inside the ComputeClaim function)
//#include <math.h>
//inline double get_claim(double r2_w2, double e, double h, double gap) { /* compute the claim */
//	double shape;
//	if (Sim.CROWNSHAPE == 0) {
//		/* CASE 0: h = h_max * (1 - (r / r_max)^exp) */
//		// BUG HERE ? should be sqrt(r2_w2), or the formula in the above comment is wrong
//		shape = 1.0 - exp(log(r2_w2) * e / 2.0);
//	} else {
//       /* CASE 1: h = h_max * (1 - r / r_max)^(1 / exp) */
//		if (e>=100) {
//			shape = 1.0;
//		} else {
//			shape = exp(log(1.0 - sqrt(r2_w2)) / e);
//		}
//	}
//	return h * (gap + (1.0 - gap) * shape);
//}


#endif





// // // //     int HSCALE = 250; /* height profile graph scale (dm) NO LESS THAN THE HIGHEST SPECIES */
// // // //       int CYCLES = 1500; /* simulation length */
// // // //       int XMAX = 512; /* should be 2^n */
// // // //       int YMAX = 512; /* should be 2^n */
// // // //       int ZMAX = 512; /* (in dm) should be at least HSCALE, perhaps more if there is a slope */
// // // //       int TMAX = 7000; /* maximum number of trees */
// // // //       int SECTORS = 8; /* number of sectors per tree: 0: off; 4, 8: on */
// // // //       double CORIENTATION = (0.5*M_PI); /* Slope orientation: highest point angle, in radians. (0.5*M_PI): North, (M_PI): West, (1.5*M_PI): South, (0): East */
// // // //       double CANGLE = (0.); /* Slope's angle (in radians, should be any value in [0,0.5*M_PI[) - 0 indicates a flat ground */
// // // //       double CSUNALT = (0.5*M_PI); /* Sun altitude (= elevation) in radians, should be any value in ]0,0.5*M_PI] - 0 indicates a sun hidden by the horizon, 0.5*M_PI indicates zenith */
// // // //       double CSUNORIENTATION = (1.54*M_PI); /* Sun azimuth. Note that this value is NOT as in the standard horizontal coordinate system, here we measure it in radians and counter clockwise. (0.5*M_PI): (light coming from) North, (M_PI): West, (1.5*M_PI): South, (0): East */
// // // //       // for crater lake (lat=42.9388312° N, long=122.1459961° W):
// // // //       // the slope angle is approx 36 degree (max altitude of the cone is 755 feet, with approx 1000 feet radius),   so ANGLE (0.2*M_PI)
// // // //       // on 21/09/2013, at 13:00 PM, elevation=47.46° (slightly more than Pi/4), azimuth=179.45° (South)
// // // //       // => zenith is for SUNORIENTATION (1.54*M_PI)   and   SUNALT (0.26*M_PI)
// // // // 
// // // //       /* ============= */
// // // //       /* MODEL OPTIONS */
// // // //       /* ============= */
// // // // 
// // // //       int SMAX = 2; /* number of species in Model */
// // // // 
// // // //       //int ROOTMAX 4 /* maximum number of trees that can have roots claim a spot */
// // // // 
// // // //       // RESTRICTCROWN
// // // //       // 0 := Off, the crown base height (tree[t].cbhs[i]) has no influence on the shape
// // // //       // 1 := the shape starts at tree[t].cbhs[i]. This will lead to a tree with a completely flat arborescence if the competition is enough to force the crown base height to "climb" to the top.
// // // //       // 2 := the shape is computed from the ground, but is cropped at tree[t].cbhs[i].
// // // //       int RESTRICTCROWN = 0; 
// // // //       // !! RESTRICTCROWN == 2 is not performing well anymore (cbhs is initialized at 2.00 in the beginning of grow_up())
// // // // 
// // // //       // VERTICALGROWTH
// // // //       // 0 := vertical growth is dependent on trunk diameter
// // // //       // 1 := vertical growth is dependent on the vertical proportion of the sector that is sunlit
// // // //       int VERTICALGROWTH = 0;
// // // // 
// // // //       // RESTRICTHEIGHT
// // // //       // this constant is used only when VERTICALGROWTH != 0
// // // //       // 0 := height increment is constant (no restriction, so tree[t].h[s] can be more than species[i].hMax)
// // // //       // 1 := height increment is constant until the tree height reaches species[i].hMax (i.e. capped)
// // // //       // 2 := height increment is proportional to | species[i].hMax - current height |
// // // //       int RESTRICTHEIGHT = 2;
// // // // 
// // // //       // RESTRICTRADIUS
// // // //       // 0 := radius increment is constant (no restriction, so tree[t].ws[s] can be more than species[i].wMax)
// // // //       // 1 := radius increment is constant until the sector radius reaches species[i].wMax (i.e. capped)
// // // //       // 2 := radius increment is proportional to | species[i].wMax - current radius at crown base |
// // // //       // 3 := radius increment is proportional to | maximal width allowed at crown base height - current radius at crown base |
// // // //       //     note that 3 is equivalent to 2 when RESTRICTCROWN == 1
// // // //       int RESTRICTRADIUS = 1;
// // // // 
// // // //       // ZRnd
// // // //       // double value allowing a tree to claim with 50% chance a tree which is ZRnd dm higher than himself
// // // //       // set at 0.0 to deactivate
// // // //       int ZRnd = 1.0;
// // // // 
// // // //       int SEEDANYWHERE = 1; /* 0 := Use newLight & newShade parameters, 1 := Seed anywhere */
// // // //       int CROWNSHAPE = 1; /* 0 := h = h_max * (1 - r^exp / r_max^exp)[Drew's profile if exp == 2]; 1 := h = h_max * (1 - r / r_max)^(1/exp)[SORTIE profile if exp == 2] */
// // // //       //int VORONOI = 0; /* 0 := Off, 1 := On */
// // // //       //int REAL3D 0 /* 0 := Off, 1 := On */
// // // // 
// // // //       /* ================ */
// // // //       /* BITMAP CONSTANTS */
// // // //       /* ================ */
// // // // 
// // // //       int OUTPUTBMP = 1; /* 0 := False, 1 := True */
// // // //       int SHOWTREELIMITS = 0; /* 0 := Hide, 1 := Show */
// // // //       //int SHOWTREEROOT = 1 /* 0 := Hide, 1 := Show */
// // // //       //int SHOWTREECENTER 1 /* 0 := Hide, 1 := Show */
// // // // //#ifdef __USE_ROOTS
// // // // //      int SHOWWATER 1 /* 0 := Hide, 1:= Show */
// // // // //#else
// // // // //      int SHOWWATER 0
// // // // //#endif


