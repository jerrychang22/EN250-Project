#ifndef _MOD_HPP_
#define _MOD_HPP_

// the following are used FOR SURE.
#include <memory>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "boost/program_options.hpp" 
#include "pngwriter.h"
#include <algorithm>
#include <iterator>
#include <unistd.h>

// the following could maybe be removed
#include <stdio.h>
#include <stdlib.h>
//#include <math.h> // conflicts with <cmath>
#include <time.h>
#include <signal.h>

// to create and remove directories:
#include <sys/stat.h>
//#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <ftw.h>
#include <unistd.h>


//#include "constants.hpp"
//#include "structs.hpp"

/* ================ */
/* BASIC STATISTICS */
/* ================ */
struct BASIC {
    double ca;                  /* sum of canopy area */
    double baCan;               /* sum of basal area of canopy trees */
    double baAll;               /* sum of basal area of all trees */
    unsigned int tNow;        /* current number of overstory trees */
    unsigned int tActive;     /* number of trees with a canopy */
    unsigned int tDead;       /* number of trees died last cycle */
    unsigned int tNew;        /* number of newborn trees  */
    unsigned int ustNow;      /* current number of understory trees */
};

struct RECORD;
struct DIB;

/* #include "structs.h" */

inline double max_ws_overstory(struct UNIT); 
inline double max_ws(struct UNIT); 
inline double min_ws(struct UNIT);
inline double get_altitude(double, double);
void initialize(std::string);
void initialize_pointers(void);
void initialize_directory(std::string);
void freeall(void);
int get_sector(int, int);
void define_header(void); /* mod_defines.cpp */
void define_species(void); /* mod_defines.cpp */
void define_colors(void); /* mod_defines.cpp */
void define_biomass(void); /* mod_defines.cpp */
void define_sectors(void); /* mod_defines.cpp */
void plant_trees(void); /* mod_growth.cpp */
void ppa(void); /* mod.cpp */
void grow_up(void); /* mod_growth.cpp */
void grow_up_ustory(void); /* mod_growth.cpp */
void water(void); /* mod_growth.cpp */
void root_competition(int *, int *); /* mod_growth.cpp */
void rootmid_competition(void); /* mod_growth.cpp */
void rootbot_competition(void); /* mod_growth.cpp */
void sunlight(void); /* mod_growth.cpp */
void implicit_root_claim(void); /* mod_growth.cpp */
void circular_root_claim(void); /* mod_growth.cpp */
void measure(void); /* mod_growth.cpp */
void color_bitmap_mainroot(void); /* mod_pictures.cpp */
void color_bitmap_cnp(void); /* mod_pictures.cpp */
void color_bitmap_hgt(void); /* mod_pictures.cpp */
void color_bitmap_ustory(void); /* mod_pictures.cpp */
void color_bitmap_ustory_hgt(void); /* mod_pictures.cpp */
void color_bitmap_rootstop(void); /* mod_pictures.cpp */
void color_bitmap_rootsbot(void); /* mod_pictures.cpp */
void color_bitmap_rootsmid(void); /* mod_pictures.cpp */
void save_bmp(std::string , int); /* mod_pictures.cpp */
void save_png_crowns(int); /* mod_pictures.cpp */
void save_png(std::string , int); /* mod_pictures.cpp */
void save_data(std::string , int); /* mod_pictures.cpp */
void save_3d(int);
void save_tree(std::string , int); /* mod_pictures.cpp */
void save_alltrees(int);
void save_stats(int);
void save_space(std::string , int); /* mod_pictures.cpp */
// early attempt, to be deleted...
//void export_space(void);
void save_balance(std::string , int); /* mod_pictures.cpp */
#if VORONOI
void voronoi(int); /* mod_pictures.cpp */
#endif
void cancel_key(int);
void plant_circles(void);

void grow_up_3d(void); /* mod_growth.cpp */
void set_dxdy(short int*, short int*);


int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf);


#endif
/* vim: set ts=4 sw=4: */
