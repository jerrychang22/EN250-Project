/* vim: set ts=4 sw=4: */

//#define __USE_ROOTS

#include "mod.hpp"
#include "structs.hpp"
Simulation* Sim;

double zslope;
double a;

double ORIENTATION,ANGLE,SUNORIENTATION,SUNALT;


/*char *file_path = (char *) "D:/Work/Models/Forest/";*/    /* data folder */
/*char *file_path = (char *) "Work/Models/Forest";*/ /* Linux data folder */
std::string file_path;
std::string file_path_root =  "../data"; /* Linux data folder */
//std::string file_path_root =  "/tmp/forestdata"; /* Linux data folder */
std::string bmp_file_mainroot =  "Bmp/MainRoot_";    /* tree canopy images */
std::string bmp_file_cnp =  "Bmp/Canopy_";    /* tree canopy images */
std::string bmp_file_hgt =  "Bmp/Height_";    /* tree height profile */
std::string bmp_file_ustorycnp =  "Bmp/UnderstoryCanopy_";    /* tree canopy images */
std::string bmp_file_ustoryhgt =  "Bmp/UnderstoryHeight_";    /* tree height profile */
std::string bmp_file_rootstop =  "Bmp/RootsTop_";    /* tree roots */
std::string bmp_file_rootsmid =  "Bmp/RootsMid_";    /* tree roots */
std::string bmp_file_rootsbot =  "Bmp/RootsBot_";    /* tree roots */
std::string txt_file_hgt =  "Data/Stat_"; /* tree height for all points on grid */
std::string txt_file_space =  "Data/Space_";    /* some tree parameters  */
std::string txt_file_tree =  "Data/Tree_";    /* some tree parameters  */
std::string dat_file_tree =  "Data/Set_"; /* complete tree data for future analysis */
std::string txt_file_bal =  "Data/Balance";   /* dead and newborn trees balance */

std::string txt_file_alltrees =  "Data/Alltrees";    /* trees from all years  */
std::ofstream file_all_trees; 

std::string txt_file_stats =  "Data/Stats";    /* trees from all years  */
std::ofstream file_stats; 

std::string heightarea_file =  "Data/Heigarea_";  /* data height, dbh and area of alive trees and their areas for all points on grid */

/* Biomasses & Profiles */
int **biomass; /* main DATA array */
unsigned int **profile; /* main HEIGHT array (mm) */
int **ustory_biomass; /* understory DATA array */
unsigned int **ustory_profile; /* understory HEIGHT array (mm) */

#ifdef __USE_ROOTS
/* Root layers and associated values */
unsigned int ***roots_top; /* highest level of roots */
const double top_water = 2.5;
unsigned int ***roots_mid; /* middle level of roots */
const double mid_water = 1.25;
const double midqualifier = 5.0; /* dbh for trees to reach mid layer */
unsigned int ***roots_bot; /* bottom level of roots */
const double bot_water = 0.5;
const double botqualifier = 10.0; /* dbh for trees to reach bot layer */
#endif


double hOver = 0;             /* projected crown cross-section overlap height */
unsigned int tTotal = 0;      /* total number of trees planted (unique tree index) */

/* ============= */
/* MODEL OPTIONS */
/* ============= */
//const unsigned int Sim->nSeedTotal = 50 + Sim->SMAX - 1;    /* total number of seed attempts per cycle (see lines 534-535) */

const double hNewRnd = 0.5;     /* newborn tree height variety 50% */
double hRnd = 0.0;        /* crown height profile variety 0% */

struct PIXEL *color;

std::vector< std::vector<struct PIXEL> > bitmap; /* main BITMAP array */
std::vector< std::vector<struct PIXEL> > bitmap_root; /* root BITMAP array */

/* set of preset colors */
const struct PIXEL pBlack = { 0x00, 0x00, 0x00 };
const struct PIXEL pWhite = { 0xFF, 0xFF, 0xFF };
const struct PIXEL pLightGrey = {0x99, 0x99, 0x99}; /* Mainly used with roots */
const struct PIXEL pGrey = {0x66, 0x66, 0x66}; /* Mainly used with roots */
const struct PIXEL pDarkGrey = {0x33, 0x33, 0x33}; /* Mainly used with roots */
const struct PIXEL pRed = { 0x00, 0x00, 0xFF };
const struct PIXEL pGreen = { 0x00, 0xFF, 0x00 };
const struct PIXEL pBlue = { 0xFF, 0x00, 0x00 };
const struct PIXEL pYellow = { 0x00, 0xFF, 0xFF };
const struct PIXEL pCyan = { 0xFF, 0xFF, 0x00 };
const struct PIXEL pPurple = { 0xFF, 0x00, 0xFF };
const struct PIXEL pLeaves = { 0x00, 0xFF, 0xC0 }; 

/*=================*/
/* TREE PARAMETERS */
/*=================*/
// struct UNIT *tree;
//std::vector<struct UNIT *> tree(Sim->TMAX);
//std::vector<struct UNIT *> tree;
//std::vector<std::unique_ptr<UNIT>> tree(Sim->TMAX);
std::vector<std::unique_ptr<Tree>> tree;
//tree.reserve(Sim->TMAX);

#ifdef __USE_ROOTS
struct GROUND *roots;
#endif

/* ================= */
/* SPECIES CONSTANTS */
/* ================= */
//struct Species *species;
inline double height(Tree tree) {   /* tree height */
    return tree.species->hMax * (1.0 - (exp(-100 * tree.d / tree.species->hMax)));
}

/* ================= */
/*        MAIN       */
/* ================= */
int main(int argc, char ** argv) {
    int i = 0;
#ifdef __USE_ROOTS
    int mid = 0, bot = 0;
#endif
    Sim = new Simulation(argc, argv);
    bitmap.resize(Sim->YMAX, std::vector<struct PIXEL>(Sim->XMAX));
    bitmap_root.resize(Sim->YMAX, std::vector<struct PIXEL>(Sim->XMAX));
    ORIENTATION = Sim->CORIENTATION;
    ANGLE = Sim->CANGLE;
    SUNORIENTATION = Sim->CSUNORIENTATION;
    SUNALT = Sim->CSUNALT;

    //std::string output_dirname;
    //if (argc>1) {
    //    output_dirname = argv[1];
    //} else {
    //    output_dirname = "default";
    //}
    initialize(Sim->output_dirname);
    //initialize((argc>1 ? argv[1] : "default"), (int) time(NULL));

    if (Sim->OUTPUTALLTREES) {
        std::string name = file_path + txt_file_alltrees + ".csv";
        file_all_trees.open(name);
        std::cout << "Just opened: " << name << "\n";
        file_all_trees << "YEAR\tID\tX\tY\tSPECIES\tDBH\tHEIGHT\tCA\tWATER_RATIO\tX0\tY0\tCBHS_1\tCBHS_2\tCBHS_3\tCBHS_4\tCBHS_5\tCBHS_6\tCBHS_7\tCBHS_8\tWS_1\tWS_2\tWS_3\tWS_4\tWS_5\tWS_6\tWS_7\tWS_8\t" << std::endl;
        file_all_trees << std::setprecision(10);
    }

    if (Sim->OUTPUTSTATS) {
        std::string name = file_path + txt_file_stats + ".csv";
        file_stats.open(name);
        std::cout << "Just opened: " << name << "\n";
        file_stats << "YEAR\tPOP\tAGE\tBA\tSTI_NB\tSTI_BA\tASSYM\tCANOPY_RATIO\tUNDER_NB\tOVER_NB\tUNDER_BA\tOVER_BA\tUNDER_AGE\tOVER_AGE\tWATER_CONSUMED\tWATER_NEEDED\tWATER_CONSUMED_A\tWATER_CONSUMED_B\tWATER_CONSUMED_C\tWATER_CONSUMED_D\tWATER_CONSUMED_E\tWATER_NEEDED_A\tWATER_NEEDED_B\tWATER_NEEDED_C\tWATER_NEEDED_D\tWATER_NEEDED_E\tWATER_CONSUMED_INTOLERANT\tWATER_NEEDED_INTOLERANT\tWATER_CONSUMED_TOLERANT\tWATER_NEEDED_TOLERANT\tPOP_A\tPOP_B\tPOP_C\tPOP_D\tPOP_E\tPOP_TOLERANT\tPOP_INTOLERANT" << std::endl;
        file_stats << std::setprecision(10);
    }

    for (i = 1; i <= Sim->CYCLES; i++) {

        //ORIENTATION = 2.0*M_PI*((double)i/(1.0+Sim->CYCLES));
        //SUNORIENTATION = M_PI * ( 2.0 - (((double)Sim->CYCLES-(double)i+1.0)/((double)Sim->CYCLES)));
        //SUNORIENTATION = M_PI * ( 2.0 - (((double)i-1.0)/((double)Sim->CYCLES)));
        //printf("current sun orientation is %.1f Pi\n",SUNORIENTATION/M_PI);
        //double pos = ((double)i/(((double)Sim->CYCLES+1.0)/2.0));
        //if ( pos < 1.0)
        //	SUNALT = M_PI * (0.025 + 0.45 * pos);
        //else if ( pos < 2.0)
        //	SUNALT = M_PI * (0.025 + 0.45 * (2.0 - pos));
        //else
        //	printf("ERROR that should not happen\n");
        a = tan(SUNORIENTATION);
        zslope = tan(SUNALT);
        //plant_circles();

        //if (i == 1) {
        //    int s = 0;
        //tree.push_back(std::unique_ptr<Tree> (new Tree(Sim)));
        //tree.back()->species = Sim->species[s];
        //tree.back()->s = s;
        //tree.back()->n = 1;
        //tree.back()->xo = 100;
        //tree.back()->yo = 100;
        //tree.back()->x = 100;
        //tree.back()->y = 100;
        //tree.back()->d = tree.back()->species->dInit * (1 + hNewRnd * ((double) (rand() % 100) / 50 - 1));
        //tree.back()->h = height(*(tree.back()));
        //tree.back()->cbh = 0;
        //tree.back()->is_ustory = 0;
        //tree.back()->pot_radius = 2.0;
        //for (int j = 0; j < Sim->SECTORS; j++) {
        //                tree.back()->max_sectors[j] = tree.back()->ws[j] = 2.0;
        //                tree.back()->ustory_ws[j] = 0;
        //                tree.back()->sector_story[j] = 0;
        //                tree.back()->cbhs[j] = tree.back()->cbh;
        //}
        //tree.back()->tr = tree.back()->species->sunShade;
        //tree.back()->ca = 0;
        //tree.back()->brdl = 0;
        //tree.back()->brdh = 0;
        //tree.back()->age = 1;
        //}

        //std::cout << "Year " << i << "\n";
        plant_trees(); // omitted for debug
        /*	ppa(); */
        grow_up();
        grow_up_ustory();
        if (Sim->ROOT == 1) {
            circular_root_claim();
        } else if (Sim->ROOT == 2) {
            implicit_root_claim();
        }
        //measure();
//#ifdef __USE_ROOTS
//        root_competition(&mid, &bot);
//        if (mid > 0) /* If mid != 0, compete on middle */
//            rootmid_competition();
//        if (bot > 0) /* Same for bot */
//            rootbot_competition();
//        water();
//#endif
        if ((Sim->OUTPUTBMP) && (Sim->OUTPUTMOD > 0) && (i % Sim->OUTPUTMOD == 0)) {
            //color_bitmap_mainroot();
            //save_bmp(bmp_file_mainroot, i);

            //color_bitmap_cnp();
            //save_png(bmp_file_cnp, i);
            
            color_bitmap_hgt();
            save_png(bmp_file_hgt, i);

            //save_png_crowns(i);
            //save_3d(i);
            //color_bitmap_ustory_hgt();
            //save_bmp(bmp_file_ustoryhgt, i);
            //color_bitmap_ustory();
            //save_bmp(bmp_file_ustorycnp, i);
        }
//#ifdef __USE_ROOTS
//        color_bitmap_rootstop();
//        save_bmp(bmp_file_rootstop, i);
//        if (mid > 0) {
//            color_bitmap_rootsmid();
//            save_bmp(bmp_file_rootsmid, i);
//        }
//        if (bot > 0) {
//            color_bitmap_rootsbot();
//            save_bmp(bmp_file_rootsbot, i);
//        }
//#endif
#if VORONOI
        voronoi(i);
#endif

        sunlight();

//        if ((Sim->OUTPUTMOD > 0) && (i % Sim->OUTPUTMOD == 0)) {
//            // save_data(dat_file_tree, i);
//            // save_space(txt_file_space, i);
//            /*save_heightarea(heightarea_file, i); */
//            /*save_profile(txt_file_hgt,i); */
//            // save_balance(txt_file_bal, i);
//        }
        if (Sim->OUTPUTALLTREES) {
            save_alltrees(i); // save for statistics at each iteration
        }
        if (Sim->OUTPUTSTATS) {
            save_stats(i); // save for statistics at each iteration
        }
    }
    freeall(); /* We're done with the pointers, say goodbye */
    return 0;
}

/* ================ */
/* BASIC STATISTICS */
/* ================ */
//struct BASIC stats[Sim->SMAX + 1]; /* all by species (Sim->SMAX for total) */

/* ===================== */
/* SPACE CLAIM STRUCTURE */
/* ===================== */
//const unsigned int cMax = Sim->TMAX;   /* maximum number of claims per point */

double ** ground_altitude; // the altitude of the ground for every integer location of the stand

const int CstMaxClaim=50; // quick and dirty max possible claim number for each pixel for the root system
struct RootClaim {                  /* claim record for each biomass pixel of overstory tree  */
    int trees_i[CstMaxClaim];              /* tree numbers  */
    double trees_s[CstMaxClaim];           /* claim strength  */
    int claim_nb; // current claim counter
    double norm_factor; // normalization factor
};

struct claimRECORD {                 /* claim record for each biomass pixel of overstory tree  */
    unsigned int t;           /* tree number  */
    double v;                   /* altitude of a claim */
    double z;                   /* claim value relative to the base of the trunk */
    unsigned short sec; 	/* sector number */
    double rad; /* current radius bein compared. */
    int is_sunny; /* 1 if shadowed, 2 otherwise */
};

struct uclaimRECORD {                 /* claim record for each biomass pixel of understory tree */
    unsigned int t;           /* tree number  */
    double v;                   /* altitude of a claim */
    double z;                   /* claim value relative to the base of the trunk */
    unsigned short sec; 	/* sector number */
    double rad; /* current radius bein compared. */
};

// early attempt, to be deleted...
//struct SPACEOCCUP {
//	unsigned short t;           /* tree number  */
////	unsigned short species;     /* tree species  */
////	unsigned short light;       /* percentage of light received  */
//} **** space;

//struct XYHEIGHT {
//	unsigned short z;           /* altitude  */
//	unsigned short t;           /* tree number  */
//	unsigned short species;     /* tree species  */
//	unsigned short light;       /* percentage of light received  */
//} *** space;

double ** water_content;
struct RootClaim *** xyrc;
struct claimRECORD *** xyclaims;
struct uclaimRECORD *** xyuclaims;
struct RECORD **** xyclaimlist;

int *** spacetree;

#ifdef __USE_ROOTS
struct RCLAIM {                 /* claim record for each root pixel  */
    unsigned int t;           /* root number/position in array  */
    double v;                   /* claim value: dbh/(distance from root pos)  */
    short dx, dy;				/* distances from root, for this claim */
    unsigned short sec; 		/* sector number */
    double rad; 				/* current radius bein compared. */
    double rad2; 				/* current radius (dx^2+dy^2) bein compared. */
} *rootclaims;
#endif

/* ============================= */
/* INITIALIZATION/FREE FUNCTIONS */
/* ============================= */
void initialize(std::string xpname) {
    setbuf(stdout, NULL);
    if (Sim->SEED==0) {
        srand(time(NULL)+rand()+getpid());
    } else {
        srand(Sim->SEED);
    }
    initialize_directory(xpname); /* this step can be removed if no output is desired */
    initialize_pointers(); /* All those pointers need memory after all */
    define_header();
    //define_colors();
    //define_species();
    define_biomass();

    a = tan(SUNORIENTATION);
    zslope = tan(SUNALT);
    //printf("a is equal to %f; 1/a = %f\n",a,1/a);
}

void initialize_pointers() {
    /* Order is semi-important although not crucial. Malloc'ing all the memory
     * consecutively makes going through the arrays faster and allows for 
     * pointer arithmetic which would make the program faster, albeit harder
     * to read after not having seen it for a while. --Ian Cordasco
     */
    int i, j;
    int k;

    // CPP
    //tree.resize(Sim->TMAX);
    //tree = malloc(Sim->TMAX * sizeof(struct UNIT));

    biomass = (int **) malloc(Sim->YMAX * sizeof(int *));
    for (i = 0; i < Sim->YMAX; i++)
        biomass[i] = (int *) malloc(Sim->XMAX * sizeof(int));

    profile = (unsigned int **) malloc(Sim->YMAX * sizeof(unsigned int *));
    for (i = 0; i < Sim->YMAX; i++)
        profile[i] = (unsigned int *) malloc(Sim->XMAX * sizeof(unsigned int));

    ustory_biomass = (int **) malloc(Sim->YMAX * sizeof(int *));
    for (i = 0; i < Sim->YMAX; i++)
        ustory_biomass[i] = (int *) malloc(Sim->XMAX * sizeof(int));

    ustory_profile = (unsigned int **) malloc(Sim->YMAX * sizeof(unsigned int *));
    for (i = 0; i < Sim->YMAX; i++)
        ustory_profile[i] = (unsigned int *) malloc(Sim->XMAX * sizeof(unsigned int));


    if (Sim->ROOT == 1) {
        // we don't bother to allocate unless we are not using an explicit root model
        water_content = (double **) malloc(Sim->XMAX * sizeof(double *));
        for (i = 0; i < Sim->XMAX; i++) {
            water_content[i] = (double *) malloc(Sim->YMAX * sizeof(double));
            for (j = 0; j < Sim->YMAX; j++) {
                if (Sim->WATER_AVAILABLE >= 0) { 
                    water_content[i][j] = Sim->WATER_AVAILABLE;
                } else {
                    // old way:
                    //WA < 0 => the water content is 0 at the center, and goes up to to the -WA

                    //water_content[i][j] = Sim->WATER_AVAILABLE * ((double)abs(i-((double)Sim->XMAX)/2.0)/(double)(Sim->XMAX/2.0) * -1.0);
                    // new way:
                    // 0   -> 1.5/8: -WA
                    // 1/8 ->   4/8: 0
                    // 3/8 ->   7/8: -WA * (i-3/8)/(4/8)
                    // 6/8 ->   1  : -WA
                    double val;
                    if (((double)i) < 1.5*((double)Sim->XMAX)/8.0  ) {
                        val = - Sim->WATER_AVAILABLE;
                    } else if (((double)i) < 4.0*((double)Sim->XMAX)/8.0  ) {
                        val = 0.0;
                    } else if (((double)i) < 7.0*((double)Sim->XMAX)/8.0  ) {
                        val = - Sim->WATER_AVAILABLE * ((double) i  - 4.0*((double)Sim->XMAX)/8.0 ) / (3.0*((double)Sim->XMAX) / 8.0);
                    } else {              
                        val = - Sim->WATER_AVAILABLE;
                    }
                    if (j==0) {
                        std::cout << i << ", " << val << std::endl;
                    }
                    water_content[i][j] = val;
                }
            }
        }
    }


    xyrc = (struct RootClaim ***) malloc(Sim->XMAX * sizeof(struct RootClaim **));
    for (i = 0; i < Sim->XMAX; i++) {
        xyrc[i] = (struct RootClaim **) malloc(Sim->YMAX * sizeof(struct RootClaim *));
        for (j = 0; j < Sim->YMAX; j++) {
            xyrc[i][j] = (struct RootClaim *) malloc(sizeof(struct RootClaim));
            for (int k=0; k<CstMaxClaim; k++) {
                xyrc[i][j]->trees_i[k] = -1;
                xyrc[i][j]->claim_nb = 0;
                xyrc[i][j]->norm_factor = 0;
            }
        }
    }

    xyclaims = (struct claimRECORD ***) malloc(Sim->XMAX * sizeof(struct claimRECORD **));
    for (i = 0; i < Sim->XMAX; i++) {
        xyclaims[i] = (struct claimRECORD **) malloc(Sim->YMAX * sizeof(struct claimRECORD *));
        for (j = 0; j < Sim->YMAX; j++) {
            xyclaims[i][j] = (struct claimRECORD *) malloc(sizeof(struct claimRECORD));
        }
    }

    xyuclaims = (struct uclaimRECORD ***) malloc(Sim->XMAX * sizeof(struct uclaimRECORD **));
    for (i = 0; i < Sim->XMAX; i++) {
        xyuclaims[i] = (struct uclaimRECORD **) malloc(Sim->YMAX * sizeof(struct uclaimRECORD *));
        for (j = 0; j < Sim->YMAX; j++) {
            xyuclaims[i][j] = (struct uclaimRECORD *) malloc(sizeof(struct uclaimRECORD));
        }
    }

    //	could be needed if we want to allow several claims at several height on the same xy square
    //	xyclaimlist = (struct RECORD ****) malloc(Sim->XMAX * sizeof(struct RECORD ***));
    //	for (i = 0; i < Sim->XMAX; i++) {
    //		xyclaimlist[i] = malloc(Sim->YMAX * sizeof(struct RECORD **));
    //		for (j = 0; j < Sim->YMAX; j++) {
    //			xyclaimlist[i][j] = malloc(CLAIMMAX * sizeof(struct RECORD *));
    //			for (k = 0; k < CLAIMMAX; k++) {
    //				xyclaimlist[i][j][k] = (struct RECORD *) malloc(sizeof(struct RECORD));
    //				//xyclaimlist[i][j][k] = (struct RECORD *) calloc(1,sizeof(struct RECORD));
    //			}
    //		}
    //	}

    // early attempt to be deleted
    //	space = (struct SPACEOCCUP ****) malloc(Sim->XMAX * sizeof(struct SPACEOCCUP ***));
    //	for (i = 0; i < Sim->XMAX; i++) {
    //		printf("%d\n",i);
    //		space[i] = malloc(Sim->YMAX * sizeof(struct SPACEOCCUP **));
    //		for (j = 0; j < Sim->YMAX; j++) {
    //			printf("%d\n",j);
    //			space[i][j] = malloc(Sim->ZMAX * sizeof(struct SPACEOCCUP *));
    //			for (k = 0; k < Sim->ZMAX; k++) {
    //				space[i][j][k] = (struct SPACEOCCUP *) malloc(sizeof(struct SPACEOCCUP));
    //				//space[i][j][k] = (struct SPACEOCCUP *) malloc(1,sizeof(struct SPACEOCCUP));
    //			}
    //		}
    //	}

    if (Sim->REAL3D) {
        printf("allocating space:\n");
        spacetree = (int ***) malloc(Sim->XMAX * sizeof(int **));
        for (i = 0; i < Sim->XMAX; i++) {
            printf("%d/100 ...",100*i/Sim->XMAX);
            spacetree[i] = (int **) malloc(Sim->YMAX * sizeof(int *));
            for (j = 0; j < Sim->YMAX; j++) {
                spacetree[i][j] = (int *) malloc(Sim->ZMAX * sizeof(int));
                for (k = 0; k < Sim->ZMAX; k++) {
                    spacetree[i][j][k] = 0;
                }
            }
        }
        printf("done.\n");
    }

    ground_altitude = (double **) malloc(Sim->XMAX * sizeof(double *));
    for (i = 0; i < Sim->XMAX; i++) {
        ground_altitude[i] = (double *) malloc(Sim->YMAX * sizeof(double));
        for (j = 0; j < Sim->YMAX; j++) {
            ground_altitude[i][j] = get_altitude((double)i,(double)j);
        }
    }

}

void freeall() {
    //int i;
    std::cout << "FIXME: memory freeing";
    //#ifdef __USE_ROOTS
    //	int j;
    //	free(rootclaims);
    //	for (i = 0; i < Sim->YMAX; i++) {
    //		for (j = 0; j < Sim->XMAX; j++)
    //			free(roots_top[i][j]);
    //		free(roots_top[i]);
    //	}
    //	free(roots_top);
    //	for (i = 0; i < Sim->YMAX; i++) {
    //		for (j = 0; j < Sim->XMAX; j++)
    //			free(roots_mid[i][j]);
    //		free(roots_mid[i]);
    //	}
    //	free(roots_mid);
    //	for (i = 0; i < Sim->YMAX; i++) {
    //		for (j = 0; j < Sim->XMAX; j++)
    //			free(roots_bot[i][j]);
    //		free(roots_bot[i]);
    //	}
    //	free(roots_bot);
    //	free(roots);
    //#endif
    //	free(color);
    //	free(species);
    //	for (i = 0; i < Sim->YMAX; i++)
    //		free(ustory_profile[i]);
    //	free(ustory_profile);
    //	for (i = 0; i < Sim->YMAX; i++)
    //		free(ustory_biomass[i]);
    //	free(ustory_biomass);
    //	for (i = 0; i < Sim->YMAX; i++)
    //		free(profile[i]);
    //	free(profile);
    //	for (i = 0; i < Sim->YMAX; i++)
    //		free(biomass[i]);
    //	free(biomass);
    //	free(tree);
    //	return;
}

int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
{
    // this routine is used in the suppression of directories
    int rv = remove(fpath);
    if (rv)
        perror(fpath);
    return rv;
}

void initialize_directory(std::string xpname) {
    std::string nameBmp;
    std::string nameData;
    char input;
    bool redo;
    /* create file_path */
    file_path = file_path_root + "/" + xpname;
    nameBmp = file_path + "Bmp";
    nameData = file_path + "Data";
    /*    char buffer[3]; */
    if (mkdir(nameBmp.c_str(),0755)!=0 || mkdir(nameData.c_str(),0755)!=0) {
        std::cerr << "Error ! " << nameBmp << " or " << nameData << " already exist or cannot be created." << std::endl;
        if (Sim->REDO_XP == 1) {
            redo = true;
        } else if (Sim->REDO_XP == 2) {
            redo = false;
        } else if (Sim->REDO_XP == 3) {
            std::ifstream statread((nameData + "/Stats.csv").c_str());
            if (!statread) {
                std::cerr << "There is no Stats.csv file for this XP! I will try to erase the XP directories:" << std::endl;
                redo = true;
            } else {
                std::string buffer;
                int c = 0;
                while (std::getline(statread, buffer)) {
                    c++;
                }
                if (c <= Sim->CYCLES) {
                    std::cerr << "Uncomplete XP; there is a Stats.csv file for this XP, and it contains " << c << " lines. I will try to erase the XP directories:" << std::endl;
                    redo = true;
                } else {
                    std::cerr << "Complete XP; there is a Stats.csv file for this XP, and it contains " << c << " lines. I am stopping here:" << std::endl;
                    redo = false;
                }
            }
        } else {
            std::cerr << "Do you want to remove these directories first? [Y/n]";
            std::cin>>input;
            if (input=='Y') {
                redo = true;
            } else {
                redo = false;
            }
        }
        if (redo) {
            // rm the dirs
            nftw(nameBmp.c_str(), unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
            nftw(nameData.c_str(), unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
            if (mkdir(nameBmp.c_str(),0755)!=0 || mkdir(nameData.c_str(),0755)!=0) {
                std::cerr << "Impossible to create " << nameBmp << " or " << nameData << "\n";
                exit(1);
            }
        } else {
            std::cerr << "Exiting\n";
            exit(0);
        }
    }
    std::cout << "Output directories " << nameBmp << " and " << nameData << " successfully created.\n";
}


/* ====================== */
/* MODEL GROWTH FUNCTIONS */
/* ====================== */
inline double diameter_inc(Tree tree) { /* dbh increase per cycle */
    return 0.5 * tree.species->growth_factor * tree.species->b1 * exp(log(tree.d) * tree.species->b2) * exp(log(tree.species->b3) * tree.d);
}

inline double crown_inc(Tree tree) { /* dbh increase per cycle */
    return 0.5 * tree.species->growth_factor * tree.species->b1 * exp(log(tree.d) * tree.species->b2) * exp(log(tree.species->b3) * tree.d)* tree.species->rinc;
}



// // // REWRITE CPP DON"T KNOW WHAT THIS IS FOR
// // // inline double pre_radius(Tree tree) {   /* projected crown radius */
// // // 	double max = max_ws_overstory(tree);
// // // 	return max + species[tree.s].rinc;
// // // 	/*return species[tree.s].wMax * (tree.h / species[tree.s].hMax) * (tree.h / species[tree.s].hMax); */
// // // }

//inline double cross_section(Tree tree, double h) {  /* projected crown cross-section at a given height */
//    /*PPA should be changed accordingly to sector approach LATER */
//    double r, gap, sum = 0;
//    short sector;
//
//    if (h >= tree.h) { /*if h is above current tree height */
//        r = 0;
//        return r;
//    }
//    else {
//        gap = species[tree.s].gap;
//        for (sector = 0; sector < Sim->SECTORS; sector++) {
//            if (h <= tree.cbhs[sector])
//                h = tree.cbhs[sector];       /* a cylinder below lowest branches */
//            if (!Sim->CROWNSHAPE)  /*if CASE 0: r = r_max * (1 - (h / h_max)^(1 / exp))    */
//                r = (h / tree.h < gap) ? tree.ws[sector] : tree.ws[sector] * exp(log(1 - (h / tree.h - gap) / (1 - gap)) / species[tree.s].exp);
//            else /*else CASE 1: r = r_max * (1 - (h / h_max))^exp */
//                r = (h / tree.h < gap) ? tree.ws[sector] : tree.ws[sector] * (1 - exp(log((h / tree.h - gap) / (1 - gap)) * species[tree.s].exp));
//            sum += r * r; /* must keep track of sector radii squared */
//        }
//        sum /= Sim->SECTORS; /*each sector's area is ( M_PI * r * r )/ 8, this compensates for that */
//    }
//    return M_PI * sum; /*return the proper area */
//}
//
inline double pre_crown(Tree tree) {    /* projected canopy area */
    if (Sim->POTENTIALCROWN == 2) {
        return(tree.pot_radius * tree.pot_radius * M_PI); // IDEALIZED POTENTIAL CROWN
    }
    if (Sim->SECTORS != 0) {
        int i;
        double sector_radius = 0;
        if (Sim->POTENTIALCROWN == 0) {
            // we assume that the current crown is what is meaningful to establish a baseline for sunlight's need
            for (i = 0; i < Sim->SECTORS; i++) {
                sector_radius += tree.max_sectors[i] * tree.max_sectors[i];
            }
            return (sector_radius * (M_PI/Sim->SECTORS)); // current crown area
        } else if (Sim->POTENTIALCROWN == 1) {
            // we assume that the potential crown based on current largest sector is what is meaningful to establish a baseline for sunlight's need
            sector_radius = tree.max_sectors[0];
            for (i = 1; i < Sim->SECTORS; i++) {
                if (sector_radius < tree.max_sectors[i]) {
                    sector_radius = tree.max_sectors[i];
                }
            }
            return (sector_radius * sector_radius * M_PI); // crown area if all sectors were equal to the current largest sector
        }
    }
    std::cerr << "pre_crown not implemented!\n" << std::endl;
    return -1;
}

inline double get_altitude(double x, double y)
{
    if (ORIENTATION >= 0.0 && ORIENTATION < M_PI_2)
        return ((x)*(fabs(cos(ORIENTATION))) + (y)*(fabs(sin(ORIENTATION)))) * sin(ANGLE);
    else if (ORIENTATION >= M_PI_2 && ORIENTATION < M_PI)
        return (((double)Sim->XMAX-x)*(fabs(cos(ORIENTATION))) + (y)*(fabs(sin(ORIENTATION)))) * sin(ANGLE);
    else if (ORIENTATION >= M_PI && ORIENTATION < 3.0 * M_PI_2)
        return (((double)Sim->XMAX-x)*(fabs(cos(ORIENTATION))) + ((double)Sim->YMAX-y)*(fabs(sin(ORIENTATION)))) * sin(ANGLE);
    else if (ORIENTATION >= 3.0 * M_PI_2 && ORIENTATION < 2.0 * M_PI)
        return (((double)x)*(fabs(cos(ORIENTATION))) + ((double)Sim->YMAX-y)*(fabs(sin(ORIENTATION)))) * sin(ANGLE);
    printf("ERROR: SLOPE ORIENTATION IS NOT IN THE [0:2 Pi[ RANGE\n");
    exit(1);
    return 0;
}




/* ================= */
/* MAX/MIN FUNCTIONS */
/* ================= */
/* /- Actually is the max of the max_sectors[] */
inline double max_ws(Tree tree) {
    /* Max of the 'max_sectors' array */
    if (Sim->SECTORS != 0) {
        int i;
        double max = tree.max_sectors[0];
        for (i = 0; i < Sim->SECTORS; i++) {
            max = (tree.max_sectors[i] > max) ? tree.max_sectors[i] : max;
            max = (tree.ws[i] > max) ? tree.ws[i] : max;
        }
        //printf("da max is %f\n",max);
        return max;
    } else {
        std::cerr << "max_ws not implemented!\n" << std::endl;
        return -1;
    }
}

/* /- Max of ONLY the overstory ws[] */
inline double max_ws_overstory(Tree tree) {
    if (Sim->SECTORS != 0) {
        int i;
        double max = tree.ws[0];
        for (i = 0; i < Sim->SECTORS; i++)
            max = (tree.ws[i] > max) ? tree.ws[i] : max;
        return max;
    }
}

inline double min_ws(Tree tree) {
    if (Sim->SECTORS != 0) {
        int i;
        double min = tree.ws[0];
        for (i = 0; i < Sim->SECTORS; i++)
            min = (tree.ws[i] < min) ? tree.ws[i] : min;
        return min;
    }
}

/* ================== */
/* BMP FILE STRUCTURE */
/* ================== */
struct DIB {                    
    /* bmp file header (ORDER MATTERS - DO NOT CHANGE EVER) */
    unsigned int bfSize;
    unsigned int bfReserved;
    unsigned int bfOffBits;
    unsigned int biSize;
    unsigned int biWidth;
    unsigned int biHeight;
    unsigned short biPlanes;
    unsigned short biBitCount;
    unsigned int biCompression;
    unsigned int biSizeImage;
    unsigned int biXPelsPerMeter;
    unsigned int biYPelsPerMeter;
    unsigned int biClrUsed;
    unsigned int biClrImportant;
} header;

/* Because DIB struct order matters, this function 
 * (define_header) is extremely important
 * Sadly bfSize is dependent upon bfOffBits, otherwise it wouldn't be.
 */
void define_header() {
    header.bfReserved = 0;      /* reserved  */
    header.bfOffBits = 0x36;    /* offset from the beginning to the bitmap data  */
    header.bfSize = header.bfOffBits + 3 * Sim->XMAX * Sim->YMAX; /* file size  */
    header.biSize = 0x28;       /* size of the BITMAPINFOHEADER structure  */
    header.biWidth = Sim->XMAX;      /* width of the image in pixels  */
    header.biHeight = Sim->YMAX;     /* height of the image in pixels  */
    header.biPlanes = 1;        /* one plane of the target device  */
    header.biBitCount = 0x18;   /* 24 bits per pixel - truecolor  */
    header.biCompression = 0;   /* no compression  */
    header.biSizeImage = 0;     /* no compression  */
    header.biXPelsPerMeter = 0xB12; /* resolution 72 dpi  */
    header.biYPelsPerMeter = 0xB12; /* resolution 72 dpi  */
    header.biClrUsed = 0;       /* number of colors - use biBitCount  */
    header.biClrImportant = 0;  /* all colors are important  */
}

//void add_species(int s)
//{
//    double b1, b2, b3;
//    species[s].hMax = 250;
//    species[s].dMax = 50;
//    species[s].wMax = 75;
//    species[s].gap = 0.1;
//    species[s].exp = 1; // 1 is for conifers
//    species[s].dInit = 0.4;
//    species[s].tiltMax = 0.01;
//    /* species[s].gRate = 0.008; */
//    // for future reference, values given in Lorimer at al (2001) for overstory trees are approx 0.5% - they rise up to 1.5% for old trees (diam>50cm)
//    species[s].mrLight = 0.01;
//
//    // for future reference, values given in Lorimer at al (2001) for understory trees are approx 1.5%
//#if defined(DEFmrShadem)
//    species[s].mrShade = 0.005; // m
//#elif defined(DEFmrShaden)
//    species[s].mrShade = 0.01;  // n
//#elif defined(DEFmrShadeo)
//    species[s].mrShade = 0.02;  // o
//#else
//    species[s].mrShade = 0.02;  // default
//#endif
//
//    species[s].newLight = 1;
//    species[s].newShade = 0.01;
//    // original value: species[s].sunShade = 0.09;
//    //species[s].sunShade = 0.09;
//    species[s].sunShade = 0.045; // suffix 'slowunder'
//    species[s].nSeed = Sim->nSeedTotal / 2;
//
//#if defined(DEFrincl)
//    species[s].rinc = 0.5;  // l
//#elif defined(DEFrincm)
//    species[s].rinc = 1;  // m
//#elif defined(DEFrincn)
//    species[s].rinc = 2;  // n
//#elif defined(DEFrinco)
//    species[s].rinc = 3;  // o
//#elif defined(DEFrincp)
//    species[s].rinc = 4;  // p
//#else
//    species[s].rinc = 2;  // default
//#endif
//
//    species[s].hinc = 5; // plausible order of magnitude, see figs. 100 and 101, p. 199 of TPOFYS
//    species[s].rootinc = 1;
//    /* Very dependent */
//    species[s].waterdep[0] = 4;
//    species[s].waterdep[1] = 0.5;
//    species[s].waterdep[2] = 2;
//    /* PRI parameters for dbh in cm */
//    b1 = 3.856766;
//    b2 = -0.553002;
//    b3 = 0.979675;
//    /* PRI parameters for dbh in dm */
//    species[s].b1 = b1 * exp(log(10) * b2) / 5;
//    species[s].b2 = b2 + 1;
//    species[s].b3 = exp(log(b3) * 10);
//    /* height map leaf color GREEN */
//    species[s].leaf.R = 0xB0;
//    species[s].leaf.G = 0xFF;
//    species[s].leaf.B = 0x00;
//#if SHOWWATER
//    /* height map water dependent leaf color blue/green */
//    species[s].waterleaf.R = 0xB0;
//    species[s].waterleaf.G = 0xFF;
//    species[s].waterleaf.B = 0xCC;
//#endif
//}
//
//void define_species() {
//
//    /*============*/
//    /* SPECIES #1 */
//    /* (shade-tolerant) */
//    /*============*/
//    add_species(0);
//    //species[0].exp = 20; // somewhat flat
//    species[0].exp = 2; // round
//    //species[0].wMax = 30;
//    /* height map leaf color YELLOW */
//    species[0].leaf.R = 0xFF;
//    species[0].leaf.G = 0xFF;
//    species[0].leaf.B = 0x00;
//    species[0].mrShade=0.01; // shade tolerant
//    species[0].rinc=1; // slow growth (radius)
//    species[0].hinc=4; // slow growth (height)
//
//    if (Sim->SMAX > 1) {
//        /*============*/
//        /* SPECIES #2 */
//        /* (shade-intolerant) */
//        /*============*/
//        add_species(1);
//        //species[1].exp = 100; // really flat
//        //species[1].exp = 20; // somewhat flat
//        species[1].exp = 2; // round
//        //species[1].wMax = 30;
//        /* height map leaf color RED */
//        species[1].leaf.R = 0xFF;
//        species[1].leaf.G = 0x00;
//        species[1].leaf.B = 0x00;
//        species[1].mrShade=0.10; // shade intolerant
//        species[1].rinc=1.25; // fast growth (radius)
//        species[1].hinc=7.5; // fast growth (height)
//#if SHOWWATER
//        /* height map water dependent leaf color blue/green */
//        species[1].waterleaf.R = 0xFF;
//        species[1].waterleaf.G = 0x84;
//        species[1].waterleaf.B = 0x00;
//#endif
//    }
//
//    if (Sim->SMAX > 2) {
//        /*============*/
//        /* SPECIES #3 */
//        /*============*/
//        add_species(2);
//        /* height map leaf color RED */
//        species[2].leaf.R = 0xFF;
//        species[2].leaf.G = 0x00;
//        species[2].leaf.B = 0x00;
//        /* height map water dependent leaf color blue/green */
//#if SHOWWATER
//        species[2].waterleaf.R = 0xFF;
//        species[2].waterleaf.G = 0x84;
//        species[2].waterleaf.B = 0x00;
//#endif
//    }
//}

//void define_colors() { // STR - color should not depend on Sim->TMAX
//	/* crown colors - now random */
//	unsigned int t;
//
//	for (t = 0; t < Sim->TMAX; t++) {
//		color[t].R = rand() % 0x100;
//		color[t].G = rand() % 0x100;
//		color[t].B = rand() % 0x100;
//	}
//	/* vacant space - now comes in pale yellow */
//	color[Sim->TMAX].R = 0xff;
//	color[Sim->TMAX].G = 0xff;
//	color[Sim->TMAX].B = 0xc0;
//}

void define_biomass() {
    unsigned short x, y, s;

    for (x = 0; x < Sim->XMAX; x++)
        for (y = 0; y < Sim->YMAX; y++) {
            biomass[y][x] = -1;   /* all space is vacant */
            profile[y][x] = 0;  /* zero height everythere */
            ustory_biomass[y][x] = -1;   /* all space is vacant */
            ustory_profile[y][x] = 0;  /* zero height everythere */
            /* set all root points to -1 so we know it is not taken */
#ifdef __USE_ROOTS
            for (s = 0; s < ROOTMAX; s++) {
                roots_top[y][x][s] = Sim->TMAX;
                roots_mid[y][x][s] = Sim->TMAX;
                roots_bot[y][x][s] = Sim->TMAX;
            }
#endif
        }
    for (s = 0; s <= Sim->SMAX; s++) {
        //stats[s].tActive = stats[s].tDead = stats[s].tNew = stats[s].tNow = 0;
        //stats[s].ustNow = 0;
    }
}

/* void swap_trees(unsigned short tree_number) {
 *     struct GROUP t;
 *     t = tree[stats[Sim->SMAX].tNow + stats[Sim->SMAX].ustNow];
 *     tree[stats[Sim->SMAX].tNow + stats[Sim->SMAX].ustNow] = tree[tree_number];
 *     tree[tree_number] = t;
 *     stats[Sim->SMAX].tNow--;
 * }
 */

void plant_trees() {
    int s, t, i;
    int xPlant, yPlant;
    double overall_mr, tTry, prob_water_mr;
    /* unsigned short dead_trees[Sim->TMAX]; */
    /* unsigned short dt, shift; */

    if (Sim->COMPLETEDISTURBANCE > 0) {
        if (Sim->COMPLETEDISTURBANCE > ((double) (rand() % 1000) / 1000.0)) {
            // case where we have a complete disturbance
            for (t = tree.size(); t > 0;) {
                t--;
                tree.erase(tree.begin()+t);
            }
        }
    }

    for (t = tree.size(); t > 0;) {
        t--;
        tree[t]->age += 1; // advance all tree by one year
        if ((Sim->ROOT >= 1) && (tree[t]->water_uptake < tree[t]->water_needed)) {
            // add a mortality probability to water-limited trees if we are using roots. Old version:
            //overall_mr = tree[t]->mr + tree[t]->water_mr * (1.0 - tree[t]->water_uptake / tree[t]->water_needed);
            prob_water_mr = tree[t]->species->water_mr * (1.0 - tree[t]->water_uptake / tree[t]->water_needed);
            if (Sim->LIGHTWATERMORTALITY == 0) {
                // only canopy position-related mortality
                overall_mr = tree[t]->mr;
            } else if (Sim->LIGHTWATERMORTALITY == 1) { // additive (wrong code, not to be used - should even remove!)
                overall_mr = tree[t]->mr + prob_water_mr;
            } else if (Sim->LIGHTWATERMORTALITY == 3) { // added probability of death
                overall_mr = tree[t]->mr + prob_water_mr - (tree[t]->mr * prob_water_mr);
            } else { // LIGHTWATERMORTALITY==2 => law of the "maximum"
                if (tree[t]->mr > prob_water_mr) {
                    overall_mr = tree[t]->mr;
                } else {
                    overall_mr = prob_water_mr;
                }
            }
        } else {
            overall_mr = tree[t]->mr;
        }
//        if ((tree[t]->ca_ratio > 0.5) && (tree[t]->dbh < 10) && (tree[t]->species->shadetolerance == 1)) {
//            overall_mr += 0.1;
//        }
        if (overall_mr > ((double) (rand() % 1000) / 1000.0)) {    /* memento mori */
            /* dead_trees[dt++] = t */
            //stats[tree[t]->s].tDead++;
            //stats[Sim->SMAX].tDead++;
            if (tree[t]->is_ustory == 0) {
                //stats[tree[t]->s].tNow--;
                //stats[Sim->SMAX].tNow--;
            }
            else {
                //stats[tree[t]->s].ustNow--;
                //stats[Sim->SMAX].ustNow--;
            }
            /* swap_trees(t); */
            //std::cout << "FIXME: inefficient moving around of trees\n";
            //if (Sim->INTERMEDIATESPECIES == 1) { // Note to future self: removing items in the middle breaks absolute addressing in the deque, so if you want to free this tiny bit of memory you need to rewrite Sim->species as a std::list
            //    Sim->species.erase(Sim->species.begin()+tree[t]->s);
            //}
            tree.erase(tree.begin()+t);
            //			for (tt = t + 1; tt < num_trees; tt++) {
            //				/* shift trees & roots up to close the gap */
            //				tree[tt - 1] = tree[tt];
            //#ifdef __USE_ROOTS
            //				roots[tt - 1] = roots[tt];
            //#endif
            //			}
        }
    }
    /* for (t = 0, tt = 0, shift = 0; t < num_trees; t++) {
     *         if ((tt < dt) && (t == dead_trees[tt])) {
     *                 tt += 1;
     *                 shift += 1;
     *                 while((t + shift) == dead_trees[tt]) {
     *                         tt += 1;
     *                         shift += 1;
     *                 }
     *         }
     *         if (shift > 0)
     *                 tree[t] = tree[t + shift];
     * }
     */       

    //for (s = 0; ((s < Sim->SMAX) || ((Sim->INTERMEDIATESPECIES == 1) && s < 2)); s++) {
    for (s = 0; s < Sim->SMAX; s++) {
        tTry = 0;
        while(tTry < Sim->species[s]->nSeed) {
            /* make nSeed attempts to plant a baby tree */
            tTry++;
            if (tTry > Sim->species[s]->nSeed) {
                // this takes care of uneven number of seeds, by adding a new seed in a probabilistic fashion
                // for example, let's say that we have: nSeed== 1.3
                // the first try has tTry==1 at this point, so tTry < nSeed => we don't enter this block
                // the second try has tTry==2 at this point, so we enter this statement.
                // we then exit the loop with probability tTry-nSeed (i.e. 70%), leaving a probability of 30% of adding a new seed
                if (tTry - Sim->species[s]->nSeed > ((double) (rand() % 1000) / 1000.0)) {
                    break;
                }
            }
            /* random position */
            xPlant = rand() % Sim->XMAX;
            yPlant = rand() % Sim->YMAX;
            /* if nothing grows there OR if we plant anywhere */
            if (Sim->SEEDANYWHERE || ((biomass[yPlant][xPlant] == -1) && ((double) rand() / RAND_MAX < Sim->species[s]->newLight)) || ((double) rand() / RAND_MAX < Sim->species[s]->newShade)) {
                tTotal++;
                //stats[s].tNew++;
                //stats[s].tNow++;
                //stats[Sim->SMAX].tNew++;
                //stats[Sim->SMAX].tNow++;
                t = tree.size();  /* create a new tree */
                if (biomass[yPlant][xPlant] == -1)
                    biomass[yPlant][xPlant] = t;    /* take the spot */
                else
                    ustory_biomass[yPlant][xPlant] = t; // WTF - does this mean that the new tree automatically outperform the other trees in understory??
#ifdef __USE_ROOTS
                /* set the index of the tree and its roots */
                roots[t].n = tree[t]->n = tTotal;
                roots[t].s = tree[t]->s = s;
                /* set the coordinates */
                roots[t].x = tree[t]->xo = xPlant;
                roots[t].y = tree[t]->yo = yPlant;
                roots[t].radius /*= tree[t]->w*/ = 2.0; /*originally= preradius(tree) */
#else

                tree.push_back(std::unique_ptr<Tree> (new Tree(Sim)));
                if (Sim->INTERMEDIATESPECIES >= 1) {
                    // generate a new random species
                    if (Sim->INTERMEDIATESPECIES == 1) {
                        Sim->species.push_back(std::shared_ptr<Species> (new Species((double) rand() / RAND_MAX)));
                    } else if (Sim->INTERMEDIATESPECIES == 2) {
                        Sim->species.push_back(std::shared_ptr<Species> (new Species((double)1.0 + (double) rand() / RAND_MAX)));
                    } else {
                        Sim->species.push_back(std::shared_ptr<Species> (new Species((double)2.0 + (double) rand() / RAND_MAX)));
                    }
                    Sim->species.back()->water_mr = Sim->species[s]->water_mr; // now correct for the water mortality and need, based on the original species (dont care which one, as it is the same in the simulations)
                    Sim->species.back()->water_needed = Sim->species[s]->water_needed;
                    tree.back()->species = Sim->species.back();
                    tree.back()->s = Sim->species.size()-1;
                } else {
                    // use the pre-defined species
                    tree.back()->species = Sim->species[s];
                    tree.back()->s = s;
                }

                tree.back()->n = tTotal;
                /* set the coordinates */
                tree.back()->xo = xPlant;
                tree.back()->yo = yPlant;
#endif
                /* no initial inclination (tilt) */
                tree.back()->x = xPlant;
                tree.back()->y = yPlant;
                /* some variety of initial sizes  */
                tree.back()->d = tree.back()->species->dInit * (1 + hNewRnd * ((double) (rand() % 100) / 50 - 1));
                /* correct height and crown radius. Since they all start off at the same height, the different sections should be the same. */

                // REWRITE CPP
                tree.back()->h = height(*(tree.back()));
                /*initial potential/super crown */
                /* crown base set at the defined height */
                /*tree.back()->cbh = (species[tree[t]->s].gap) * tree[t]->h; */
                tree.back()->cbh = 0;
                tree.back()->is_ustory = 0;
                tree.back()->pot_radius = 2.0;
                if (Sim->SECTORS != 0)
                    for (i = 0; i < Sim->SECTORS; i++) {
                        tree.back()->max_sectors[i] = tree.back()->ws[i] = 2.0;
                        tree.back()->ustory_ws[i] = 0;
#ifdef __USE_ROOTS
                        roots.back().radii_top[i] = 2.0;
                        roots.back().radii_mid[i] = 1.0; 
                        roots.back().radii_bot[i] = 1.0;
#endif
                        tree.back()->sector_story[i] = 0;
                        tree.back()->cbhs[i] = tree.back()->cbh;
                    }
                /* all grow up equally */
                /* tree.back()->gr = species[s].gRate; */
                /* all are equally transparent */
                tree.back()->tr = tree.back()->species->sunShade;
#ifdef __USE_ROOTS
                for (i = 0; i < ROOTMAX; i++)
                    tree.back()->wrealuptake[i] = 0;
                tree.back()->wlimit = 0;
#endif
                /* no crown area - will fight for it */
                tree.back()->ca = 0;
                /* zero crown border length */
                tree.back()->brdl = 0;
                /* zero crown border average height */
                tree.back()->brdh = 0;
                /* one year old */
                tree.back()->age = 1;
                /* canopy color */
                //tree.back()->color = color[t];
                tree.back()->color.R = rand() % 0x100;
                tree.back()->color.G = rand() % 0x100;
                tree.back()->color.B = rand() % 0x100;

                //// NEW ADDITION WHILE REWRITE CPP
                // Is this really an optimization??
                //tree.back()->claims = (double **) malloc(species[s].wmax * sizeof(double *));
                //for (int ii = 0; ii < species[s].wmax; ii++) {
                //    tree.back()->claims[ii] = (double *) malloc(species[s].wmax * sizeof(double));
                //    for (jj = 0; jj < species[s].wmax; jj++) {
                //        tree.back()->claims[ii][jj] = 0;
                //    }
                //}
            }
        }
        //printf("+%3i[%3i] ", stats[s].tNew, tTry);
    }
}

//void ppa() {
//    double h, hMin, hMax, area;
//    unsigned int num_trees = tree.size();
//    unsigned int t;
//    /* projected crown cross-section overlap height */
//    hMin = 0;
//    hMax = Sim->HSCALE;
//    while((hMax - hMin) > 0.01) {   /* precision of 1 mm */
//        area = 0;
//        h = (hMax + hMin) / 2;
//        for (t = 0; t < num_trees; t++)
//            // REWRITE CPP
//            area += (tree[t]->is_ustory == 0) ? cross_section(*tree[t], h) : 0;
//        if (area < Sim->XMAX * Sim->YMAX)
//            hMax = h;
//        else
//            hMin = h;
//    }
//    hOver = hMin;
//}

//early attempt, to be deleted..
//void export_space() {
//	unsigned short sector, t;
//	unsigned short num_trees = stats[Sim->SMAX].tNow + stats[Sim->SMAX].ustNow;
//	double sector_angle;
//	short x, y, z;
//	for (x = 0; x < Sim->XMAX; x++) {
//		for (y = 0; y < Sim->YMAX; y++) {
//				xyclaims[x][y]->v = 0;
//				//space[x][y][z]->species = -1;
//				//space[x][y][z]->light = 0;
//		}
//	}
//
// 	for (t = 0; t < num_trees; t++) {
//		for (sector = 0; sector < Sim->SECTORS; sector++) {
//			sector_angle = (M_PI * 2.0 * ((double) sector) + 1.0) / ((double) Sim->SECTORS);
//			x = tree[t]->x + cos(sector_angle);
//			y = tree[t]->y + sin(sector_angle);
//			if (x<0) {
//				x += Sim->XMAX;
//			} else if (x >= Sim->XMAX) {
//				x -= Sim->XMAX;
//			}
//			if (y<0) {
//				y += Sim->XMAX;
//			} else if (y >= Sim->XMAX) {
//				y -= Sim->XMAX;
//			}
//			z = tree[t]->ws[sector] * cos(sector_angle - ORIENTATION) * tan(ANGLE);
//
//			xyclaims[x][y]->z = z;
//			//space[x][y][z]->species = tree[t]->s;
//			//space[x][y][z]->light = 0;
//		}
//	}
//}


// // // DONT WANT TO BE BOTHERED BY THAT DURING CPP REWRITE
// // void grow_up_3d() {
// // 	int x, y, z, t, xi, yi;
// // 	short endx, endy, starty;
// // 	int num_trees = stats[Sim->SMAX].tNow + stats[Sim->SMAX].ustNow;
// // 	/* dx, dy MUST be signed */
// // 	short dx, dy; /* Used to determine periodic boundaries to find correct
// // 					 sector. */
// //     // CPP
// // 	//struct UNIT *temptree = (struct UNIT *)malloc(num_trees * sizeof(struct UNIT));
// //     std::vector<Tree *> temptree;
// //     for (int to_be_deleted_i=0; to_be_deleted_i<num_trees; to_be_deleted_i++) {
// //         temptree.push_back(new Tree());
// //         //temptree[to_be_deleted_i] =  (Tree*) malloc(sizeof(struct UNIT));
// //     }
// //     //temptree.reserve(num_trees);
// // 	int sector;
// // 	double sector_angle;
// // 	short delta;
// // 	//unsigned int c = 0;
// // 	//	unsigned int trees_around = 0;
// // 	double w2, r2, v;
// // 	double x2,y2,nextx,nexty;
// // 	short xn,yn;
// // 	double ax,ay;
// // 	double xsgn,ysgn;
// // 	double lighttraveldist;
// // 	double current_altx,current_alty,light_alt,xshift,yshift;
// // 
// // 	delta = 0;
// // 
// // 	for (t = 0; t < num_trees; t++) {
// // 		tree[t]->ca = 0;
// // 		if ( tree[t]->is_ustory == 2 ) {
// // 			// under and overstory tree
// // 			for (sector = 0; sector < Sim->SECTORS; sector++) {
// // 				if (tree[t]->ws[sector] > tree[t]->ustory_ws[sector]) {
// // 					tree[t]->max_sectors[sector] = tree[t]->ws[sector];
// // 				} else {
// // 					tree[t]->max_sectors[sector] = tree[t]->ustory_ws[sector];
// // 				}
// // 			}
// // 		} else {
// // 			if ( tree[t]->is_ustory == 1 ) {
// // 			// understory tree
// // 				for (sector = 0; sector < Sim->SECTORS; sector++) {
// // 					tree[t]->max_sectors[sector] = tree[t]->ustory_ws[sector];
// // 				}
// // 			} else {
// // 				/* otherwise all the largest sectors are in the overstory */
// // 				for (sector = 0; sector < Sim->SECTORS; sector++)
// // 					tree[t]->max_sectors[sector] = tree[t]->ws[sector];
// // 			}
// // 		}
// //             // REWRITE CPP
// // 		tree[t]->pot_ca = pre_crown(*tree[t]);
// // 		/*Reset sectors' radii */
// // 		for (sector = 0; sector < Sim->SECTORS; sector++) {
// // 			temptree[t]->ws[sector] = 2.00;
// // 			//temptree[t]->cbhs[sector] = tree[t]->cbhs[sector];
// // 			temptree[t]->cbhs[sector] = tree[t]->h;
// // 		}
// // 		temptree[t]->w = max_ws(*tree[t]); /* maximum value of sectors' radii */
// // 	}
// // 
// // 	stats[Sim->SMAX].tNow = stats[Sim->SMAX].ustNow = 0;
// // 
// // 	for (x = 0; x < Sim->XMAX; x++) {
// // 		for (y = 0; y < Sim->YMAX; y++) {
// // 			for (z = 0; z < Sim->ZMAX; z++) {
// // 				spacetree[x][y][z] = -1;
// // 			}
// // 		}
// // 	}
// // 
// // 	for (t = 0; t < num_trees; t++) {
// // 		for (xi = -temptree[t]->w; x <= temptree[t]->w ; x++) {
// // 			for (yi = -temptree[t]->w; y < -temptree[t]->w; y++) {
// // 
// // 				/* distance */
// // 				x = temptree[t]->x - xi;
// // 				if (x<0) {
// // 					x+=Sim->XMAX;
// // 				} else if (x >= Sim->XMAX) {
// // 					x-=Sim->XMAX;
// // 				}
// // 				y = temptree[t]->y - yi;
// // 				if (y<0) {
// // 					y+=Sim->YMAX;
// // 				} else if (y >= Sim->YMAX) {
// // 					y-=Sim->YMAX;
// // 				}
// // 				dx = x - tree[t]->x;
// // 				dy = y - tree[t]->y;
// // 
// // 				set_dxdy(&dx,&dy);
// // 
// // 				/* periodic boundaries: determines how the trees grow */
// // 				//if ((2 * abs(dx)) > Sim->XMAX) { /* if this goes over the edge of the picture */
// // 				//	if (dx > 0)
// // 				//		dx = dx - Sim->XMAX; /* dx positive */
// // 				//	/* If the point goes off the positive side of
// // 				//	 * the plot, then it becomes negative since it
// // 				//	 * is then placed to the left of the tree's
// // 				//	 * root.
// // 				//	 * Else, it is negative and becomes positive.
// // 				//	 * The same logic applies to dy.
// // 				//	 */
// // 				//	else
// // 				//		dx = Sim->XMAX + dx;
// // 				//	/* Since dx is negative, adding it to Sim->XMAX is
// // 				//	 * equivalent to Sim->XMAX - abs(dx). This is just
// // 				//	 * faster since you're not first calling 
// // 				//	 * another function.
// // 				//	 * Same applies to dy.
// // 				//	 */
// // 				//}
// // 				//if ((2 * abs(dy)) > Sim->YMAX) { /* same test as above */
// // 				//	if (dy > 0)
// // 				//		dy = dy - Sim->YMAX;
// // 				//	else
// // 				//		dy = Sim->YMAX + dy;
// // 				//}
// // 
// // 				r2 = (double) ((dx * dx) + (dy * dy));
// // 
// // 				w2 = temptree[t]->w * temptree[t]->w; /* super crown */
// // 				if (r2 < w2) { /* if point is within the projected super crown */
// // 					if (Sim->SECTORS != 0) {
// // 						sector = get_sector(dx,dy);
// // 						sector_angle = (M_PI * 2.0 * ((double) sector) + 1.0) / ((double) Sim->SECTORS);
// // 						/*Actual "crown" for comparison */
// // 						/*w2 = tree[t]->ws[sector] * tree[t]->ws[sector]; */
// // 						w2 = tree[t]->max_sectors[sector] * tree[t]->max_sectors[sector];
// // 					}
// // 					if (r2 < w2) {
// // 						/*If it is a point within that sector's "crown" */
// // 						v = get_claim(r2/w2,species[tree[t]->s].exp,tree[t]->h,species[tree[t]->s].gap) - sqrt(r2) * cos(sector_angle - Sim->CORIENTATION) * tan(Sim->CANGLE); /* the 2nd term corresponds to the change of altitude caused by the slope */
// // 						if (v<0) {
// // 							v = 0.001;
// // 						}
// // 
// // 						/* Logic of below: If you want to restrict how the crown
// // 						 * grows, then if RESTRICTCROWN = 1, !1 = 0, which forces
// // 						 * the if statement to use the v >= tree[t]->cbhs[sector]
// // 						 * If RESTRICTCROWN = 0, !0 is always true and will skip
// // 						 * the second statement.
// // 						 * -Ian Cordasco (explained)
// // 						 */
// // 						if ( !RESTRICTCROWN || v >= tree[t]->cbhs[sector] ) {
// // 							int z = v;
// // 							spacetree[x][y][z] = t;
// // 						} // end if restrict crown
// // 					}
// // 				}
// // 			}
// // 		}
// // 	}
// // 
// // 	if (SUNORIENTATION >= M_PI_2 && SUNORIENTATION < 3.0 * M_PI_2) {
// // 		xsgn = 1;
// // 		x = 0;
// // 		endx = Sim->XMAX;
// // 	} else {
// // 		xsgn = -1;
// // 		x = Sim->XMAX-1;
// // 		endx = 1;
// // 	}
// // 	if (SUNORIENTATION >= 0. && SUNORIENTATION < M_PI) {
// // 		ysgn=-1;
// // 		starty = Sim->YMAX-1;
// // 		endy = 1;
// // 	} else {
// // 		ysgn=1;
// // 		starty = 0;
// // 		endy = Sim->YMAX;
// // 	}
// // 
// // 	a = tan(SUNORIENTATION);
// // 
// // 	ax = fabs(tan(SUNORIENTATION - M_PI_2)); // slope coefficient, used for the computation of the next x value on the intersection with one integer y value
// // 	ay = fabs(a); // same for the next y value
// // 	xshift = ground_altitude[Sim->XMAX-1][Sim->YMAX-1] - ground_altitude[0][Sim->YMAX-1];
// // 	yshift = ground_altitude[Sim->XMAX-1][Sim->YMAX-1] - ground_altitude[Sim->XMAX-1][0];
// // 
// // 	//y=starty;
// // 
// // 	//#pragma omp parallel for default(none) private(x,y,x2,y2,nextx,nexty,xn,yn,lighttraveldist,xdist,ydist) shared(xsgn,ysgn,ax,ay,xyclaims,zslope,axdist,aydist)
// // 	for (; x*xsgn < endx; x += xsgn) {
// // 		for (y=starty; y*ysgn < endy; y += ysgn) {
// // 			//step 1: determine the corner that will likely cast the biggest shadow, in the current square
// // 			xn = 0; // the number of squares we are, along x axis
// // 			yn = 0; // .. .. .. .. .. .. .. .. .. .. .. .. y axis
// // 			x2 = ((double)x) + xsgn/2.0;
// // 			nextx = x2 + xsgn*ax/2.0;
// // 			y2 = ((double)y) + ysgn/2.0;
// // 			nexty = y2 + ysgn*ay/2.0;
// // 			lighttraveldist = 0;
// // 			//xdist = 0;
// // 			//ydist = 0;
// // 			current_altx = 0;
// // 			current_alty = 0;
// // 			//D printf("i will start with %d,%d\n(reminder: ax=%f and ay=%f and abs(a)=%f)\n",x,y,ax,ay,fabs(a));
// // 			//while (x+xn > 1 && y+yn > 1 && x+xn < Sim->XMAX-1 && y+yn < Sim->YMAX-1) { // stop condition if outside the boundaries
// // 			//printf("(%d,%d,xyclaims[x][y]->z - lighttraveldist*zslope | %.2f) - ",x+xn,y+yn,(xdist*cos(ORIENTATION) + ydist*sin(ORIENTATION)) * sin(ANGLE));
// // 			//while (xyclaims[x][y]->z - lighttraveldist*zslope >= + (xdist*cos(ORIENTATION) + ydist*sin(ORIENTATION)) * sin(ANGLE) ) { // stop condition is: does the lightbeam hit the ground at xn,yn ?
// // 			//while (xyclaims[x][y]->z - lighttraveldist*zslope >= get_altitude((double)x+xn,(double)y+yn) ) { // stop condition is: does the lightbeam hit the ground at xn,yn ?
// // 			light_alt = Sim->ZMAX-1;
// // 			int tree_encountered = 0;
// // 			//if (light_alt > max_alt) {
// // 			//	light_alt -= max_alt;
// // 			//} else if if (light_alt < 0) {
// // 			//	light_alt += max_alt;
// // 			//}
// // 			while (light_alt >= get_altitude(x+xn,y+yn)-1) { // stop condition is: does the lightbeam hit the ground at xn,yn ?
// // 				////step 2: compute the coordinates of the first intersection with another square.
// // 				D printf("at (%.2f,%.2f) [nextx=%.2f,nexty=%.2f]; next intersections are (%.2f,%.2f) and (%.2f,%.2f)\n",x2,y2,nextx,nexty,nextx,y+yn+ysgn,x+xn+xsgn,nexty);
// // 				if ( ysgn * (y+yn+ysgn - nexty) < 0 ) {
// // 					D printf("intersected y=%f !\n",y+yn+ysgn);
// // 					lighttraveldist +=sqrt((x2-nextx)*(x2-nextx)+(y+yn+ysgn-y2)*(y+yn+ysgn-y2));
// // 					x2 = nextx;
// // 					nextx += ax*xsgn;
// // 					//nextx = x2+ax*xsgn;
// // 					y2 = y+yn+ysgn;
// // 					yn += ysgn;
// // 					//ydist += ysgn;
// // 				} else {
// // 					D printf("intersected x=%f !\n",x+xn+xsgn);
// // 					lighttraveldist +=sqrt((y2-nexty)*(y2-nexty)+(x+xn+xsgn-x2)*(x+xn+xsgn-x2));
// // 					y2 = nexty;
// // 					nexty += ay*ysgn;
// // 					//nexty = y2+ay*ysgn;
// // 					x2 = x+xn+xsgn; //wtf
// // 					xn += xsgn;
// // 					//xdist += xsgn; 
// // 				}
// // 				//lighttraveldist = sqrt(xdist*xdist+ydist*ydist);
// // 				if (y+yn < 0) {
// // 					yn+=Sim->YMAX;
// // 					y2+=(double)Sim->YMAX;
// // 					nexty+=(double)Sim->YMAX;
// // 					current_alty--;
// // 				} else if (y+yn >= Sim->YMAX) {
// // 					yn-=Sim->YMAX;
// // 					y2-=(double)Sim->YMAX;
// // 					nexty-=(double)Sim->YMAX;
// // 					current_alty++;
// // 				}
// // 				if (x+xn < 0) {
// // 					xn+=Sim->XMAX;
// // 					x2+=(double)Sim->XMAX;
// // 					nextx+=(double)Sim->XMAX;
// // 					current_altx--;
// // 				} else if (x+xn >= Sim->XMAX) {
// // 					xn-=Sim->XMAX;
// // 					x2-=(double)Sim->XMAX;
// // 					nextx-=(double)Sim->XMAX;
// // 					current_altx++;
// // 				}
// // 				//current_alt = get_altitude((double)x+xn,(double)y+yn);
// // 				light_alt -= - lighttraveldist*zslope - current_altx * xshift - current_alty * yshift;
// // 				// step 3: delete the claim if it is below the shadowed z
// // 				//if ((xyclaims[x+xn][y+yn]->v > 0) && (xyclaims[x+xn][y+yn]->z < xyclaims[x][y]->z - lighttraveldist*zslope)) {
// // 				//if ((xyclaims[x+xn][y+yn]->z < xyclaims[x][y]->z - lighttraveldist*zslope)) {
// // 				if (spacetree[x+xn][y+yn][(int)light_alt] != -1) {
// // 					//D printf(" - just shadowed %d,%d (z=%f) - ",x-xn,y-yn,xyclaims[x+xn][y+yn]->z,x,y,xyclaims[x][y]->z);
// // 					//if (xyclaims[x][y]->is_sunny == 1) {
// // 					//						if (xyclaims[x+xn][y+yn]->is_sunny > 0) {
// // 					//							//								break; // no need to continue to look if we land on an already shaded spot
// // 					//						} else {
// // 					//							xyclaims[x+xn][y+yn]->is_sunny += 1;
// // 					//						}
// // 					if (tree_encountered) {
// // 						spacetree[x+xn][y+yn][(int)light_alt] = -2;
// // 					} else {
// // 						tree_encountered = 1;
// // 					}
// // 					//break; // no need to continue to look if the ray hit an obstacle !
// // 				//} else {
// // 				} // end if shadow claim
// // 			} // end while
// // 		} // end y loop
// // 	} // end x loop
// // 
// // 	//#pragma omp barrier
// // 
// // 	//#pragma omp parallel for default(none) private(x,y,t)
// // 	for (x = 0; x < Sim->XMAX; x++) {
// // 		for (y = 0; y < Sim->YMAX; y++) {
// // 			biomass[y][x] = -1;
// // 			profile[y][x] = 0;
// // 			for (z = Sim->ZMAX; z >= 0; z--) {
// // 				if (spacetree[x][y][z] != -1) {
// // 					t = spacetree[x][y][z]; /* the highest unshadowed tree  */
// // 					biomass[y][x] = t;
// // 					//profile[y][x] = (unsigned short) (xyclaims[x][y]->v * 100);
// // 					profile[y][x] = (unsigned int) (xyclaims[x][y]->v * 100);
// // 					tree[t]->ca++;
// // 					/*test to see if claim is greater than current cbh for sector */
// // 					/*if so, replace. */
// // 					if (Sim->SECTORS != 0) {
// // 						/* Below code checks EVERY POINT.
// // 						 * If current point is farthest point in sector
// // 						 * then you want to save it. The final longest radius
// // 						 * is then found for the tree's sector's new radius.
// // 						 */
// // 						dx = x - tree[t]->x;
// // 						dy = y - tree[t]->y;
// // 						set_dxdy(&dx,&dy);
// // 						sector = xyclaims[x][y]->sec;
// // 						if ( temptree[t]->cbhs[sector] > z-ground_altitude[x][y] ) {
// // 							temptree[t]->cbhs[sector] = z-ground_altitude[x][y];
// // 						}
// // 						if (temptree[t]->ws[sector] < xyclaims[x][y]->rad) { /* valid since an initialized sector is in squared units */
// // 							temptree[t]->ws[sector] = xyclaims[x][y]->rad;
// // 						}
// // 					} else {
// // 						tree[t]->cbh = z-ground_altitude[x][y];
// // 					}
// // 				break;
// // 				}
// // 			}
// // 		} // end y loop
// // 	} // end x loop
// // 
// // 	for (t = 0; t < num_trees; t++) {
// // 		for (sector = 0; sector < Sim->SECTORS; sector++) {
// // 			tree[t]->ws[sector] = sqrt(temptree[t]->ws[sector]); /* corrects for ws being in squared units */
// // 			tree[t]->cbhs[sector] = temptree[t]->cbhs[sector];
// // 		}
// // 	}
// // 
// // 	/*now that you have realized ws, separate overstory and understory trees */
// // 	for (t = 0; t < num_trees; t++) {
// // 		if (tree[t]->ca <= 0) {
// // 			tree[t]->is_ustory = 1; /*understory */
// // 			for (sector = 0; sector < Sim->SECTORS; sector++)
// // 				tree[t]->sector_story[sector] = 1;
// // 		}
// // 		else {
// // 			for (sector = 0; sector < Sim->SECTORS; sector++) {
// // 				if (tree[t]->ws[sector] < (tree[t]->max_sectors[sector] * 0.1)) {
// // 					/* this test aims to preserve a bit of an overstory crown
// // 					 * which just got critically overtopped (loss > 90%) */
// // 					tree[t]->sector_story[sector] = 2; /* under and overstory */
// // 				}
// // 				else {
// // 					tree[t]->sector_story[sector] = 0; /*overstory */
// // 					/*since this has reached here, kill the understory */
// // 					tree[t]->ustory_ws[sector] = 0;
// // 					tree[t]->ustory_ca = 0;
// // 				}
// // 			}
// // 			tree[t]->is_ustory = 0;
// // 			delta = 0;
// // 			for (sector = 0; sector < Sim->SECTORS; sector++)
// // 				delta += tree[t]->sector_story[sector];
// // 			/*if the sum of the sectors is 0 then the tree is purely overstory */
// // 			/*else it exists in both. */
// // 			tree[t]->is_ustory = ((delta == 0) ? 0 : 2);
// // 		}
// // 		if (tree[t]->is_ustory == 0)
// // 			stats[Sim->SMAX].tNow++;
// // 		else
// // 			stats[Sim->SMAX].ustNow++;
// // 
// // 	}
// //     // REWRITE CPP
// //     std::cout << "FIXME: memory\n";
// // 	//free(temptree);
// // }



void grow_up() {
    int x, y, t;
    short endx, endy, starty;
    int num_trees = tree.size();
    /* dx, dy MUST be signed */
    short dx, dy; /* Used to determine periodic boundaries to find correct
                     sector. */
    int sector;
    double sector_angle;
    short delta;
    //unsigned int c = 0;
    //	unsigned int trees_around = 0;
    double w2, w2_previous, r2, v, z;
    double x2,y2,nextx,nexty;
    short xn,yn;
    double ax,ay;
    double xsgn,ysgn;
    double lighttraveldist;
    double current_altx,current_alty,light_alt,xshift,yshift;

    delta = 0;

    // initialize tree[t]->max_sectors as the max width in over and understory
    // REWRITE CPP: SHOULD SIMPLIFY AND PUT THIS IN CLASS
    for (t = 0; t < num_trees; t++) {
        tree[t]->ca = 0;
        if ( tree[t]->is_ustory == 2 ) {
            // under and overstory tree
            for (sector = 0; sector < Sim->SECTORS; sector++) {
                if (tree[t]->ws[sector] > tree[t]->ustory_ws[sector]) {
                    tree[t]->max_sectors[sector] = tree[t]->ws[sector];
                } else {
                    tree[t]->max_sectors[sector] = tree[t]->ustory_ws[sector];
                }
            }
        } else {
            if ( tree[t]->is_ustory == 1 ) {
                // understory tree
                for (sector = 0; sector < Sim->SECTORS; sector++) {
                    tree[t]->max_sectors[sector] = tree[t]->ustory_ws[sector];
                }
            } else {
                /* otherwise all the largest sectors are in the overstory */
                for (sector = 0; sector < Sim->SECTORS; sector++)
                    tree[t]->max_sectors[sector] = tree[t]->ws[sector];
            }
        }
        tree[t]->pot_ca = pre_crown(*tree[t]);
        /*Reset sectors' radii */
        tree[t]->max_ws = max_ws(*tree[t]); /* maximum value of sectors' radii */
        for (sector = 0; sector < Sim->SECTORS; sector++) {
            tree[t]->ws[sector] = 2.00; // we dwarf all sectors (but they will regrow according to light claims)
            //temptree[t]->cbhs[sector] = tree[t]->cbhs[sector];
            tree[t]->cbhs[sector] = tree[t]->h; // we set the crown base height at the tree height value (but it will shrink according to light claims
        }
        //temptree[t]->w = max_ws(*tree[t]); /* maximum value of sectors' radii */
    }

    //stats[Sim->SMAX].tNow = stats[Sim->SMAX].ustNow = 0;


    //// BRANCHING HERE TO DO THE NEW VERSION TO UPDATE XYCLAIMS[X][Y]
    //	for (x = 0; x < Sim->XMAX; x++) {
    //		for (y = 0; y < Sim->YMAX; y++) {
    //			xyclaims[x][y]->v = 0;
    //			xyclaims[x][y]->z = ground_altitude[x][y];
    //        }
    //    }
    //	for (t = 0; t < num_trees; t++) {
    //        // here, update xyclaims[x,y]
    //    }
    //
    //// END BRANCHING HERE

    //#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma omp parallel for default(none) private(x,y,t,dx,dy,r2,w2,v,z,sector,sector_angle) shared(xyclaims,tree,temptree,species,num_trees,ground_altitude)
    for (x = 0; x < Sim->XMAX; x++) {
        for (y = 0; y < Sim->YMAX; y++) {
            xyclaims[x][y]->v = 0;
            xyclaims[x][y]->z = ground_altitude[x][y];
            //xyclaims[x][y]->is_sunny = 9999; // arbitrary high shadow
            xyclaims[x][y]->is_sunny = 2; // every square is shaded, by default

            for (t = 0; t < num_trees; t++) {
                /* distance */
                dx = x - tree[t]->x;
                dy = y - tree[t]->y;

                /* periodic boundaries: determines how the trees grow */
                if ((2 * abs(dx)) > Sim->XMAX) { /* if this goes over the edge of the picture */
                    if (dx > 0)
                        dx = dx - Sim->XMAX; /* dx positive */
                    /* If the point goes off the positive side of
                     * the plot, then it becomes negative since it
                     * is then placed to the left of the tree's
                     * root.
                     * Else, it is negative and becomes positive.
                     * The same logic applies to dy.
                     */
                    else
                        dx = Sim->XMAX + dx;
                    /* Since dx is negative, adding it to Sim->XMAX is
                     * equivalent to Sim->XMAX - abs(dx). This is just
                     * faster since you're not first calling 
                     * another function.
                     * Same applies to dy.
                     */
                }
                if ((2 * abs(dy)) > Sim->YMAX) { /* same test as above */
                    if (dy > 0)
                        dy = dy - Sim->YMAX;
                    else
                        dy = Sim->YMAX + dy;
                }

                r2 = (double) ((dx * dx) + (dy * dy));

                w2 = tree[t]->max_ws * tree[t]->max_ws; /* super crown */
                
                if (r2 < w2) { /* if point is within the projected super crown */
                    if (Sim->SECTORS != 0) {
                        sector = get_sector(dx,dy);
                        sector_angle = (M_PI * 2.0 * ((double) sector) + 1.0) / ((double) Sim->SECTORS);
                        /*Actual "crown" for comparison */
                        /*w2 = tree[t]->ws[sector] * tree[t]->ws[sector]; */
                        //w2_previous = tree[t]->max_sectors[sector] * tree[t]->max_sectors[sector]; // max_sectors conveniently stores the last width of the sector
                        w2_previous = 2.0+tree[t]->max_sectors[sector] * tree[t]->max_sectors[sector]; // max_sectors conveniently stores the last width of the sector
                    } else {
                        // handle the case where there are no sector
                        sector = 0;
                        w2_previous = w2;
                        sector_angle = 0;
                        sector = 0;
                    }
                                //if (sector == 6) {
                                //    std::cerr << r2 << "..." << w2_previous << std::endl;
                                //}

                    if (r2 < w2_previous) {
                        /*If it is a point within that sector's "crown" */
                        v = tree[t]->ComputeClaim(r2, w2_previous, sector);
                        z = v - sqrt(r2) * cos(sector_angle - Sim->CORIENTATION) * tan(Sim->CANGLE); /* the 2nd term corresponds to the change of altitude caused by the slope */
                        if (z<=0 && v > 0) {
                            // added to consider that branches in collision with the slope are actually bended on the ground, at an arbitrary small altitude
                            // the additional check of v>0 is to make sure that we do restore a canceled claim (because of RESTRICTCROWN = 2)
                            z = 0.001; 
                        }
                        z += ground_altitude[x][y];

                        // if ( z > xyclaims[x][y]->z || (z >= xyclaims[x][y]->z - Sim->ZRnd && (rand() % 2) && z > ground_altitude[x][y]) ) {
                        if ( z > xyclaims[x][y]->z || ((z == xyclaims[x][y]->z) && (rand() % 2) && (z > ground_altitude[x][y])) ) {
#pragma omp critical(claimupdate)
                            { // begin pragma critical
                                xyclaims[x][y]->rad = sqrt(r2); // the radius in dm of the claim for this tree
                                xyclaims[x][y]->t = t;
                                xyclaims[x][y]->sec = sector;
                                xyclaims[x][y]->v = v;
                                xyclaims[x][y]->z = z;
                                xyclaims[x][y]->is_sunny = 0; // we suppose each claim is not shadowed 
                                //printf("x=%d ; y=%d; z=%.2f; v=%.2f\n",x,y,z,v);
                            } // end pragma critical
                        } // end if x,y claimed
                    }
                }
            }
        }
    }

    if (SUNALT != 0.5*M_PI) {
        // If the sun is not at zenith, we need to compute custom tree shading
        printf("!!! Will now compute sun shading - currently this might be buggy so beware\n");

        if (SUNORIENTATION >= M_PI_2 && SUNORIENTATION < 3.0 * M_PI_2) {
            xsgn = 1;
            x = 0;
            endx = Sim->XMAX;
        } else {
            xsgn = -1;
            x = Sim->XMAX-1;
            endx = 1;
        }
        if (SUNORIENTATION >= 0. && SUNORIENTATION < M_PI) {
            ysgn=-1;
            starty = Sim->YMAX-1;
            endy = 1;
        } else {
            ysgn=1;
            starty = 0;
            endy = Sim->YMAX;
        }

        a = tan(SUNORIENTATION);

        ax = fabs(tan(SUNORIENTATION - M_PI_2)); // slope coefficient, used for the computation of the next x value on the intersection with one integer y value
        ay = fabs(a); // same for the next y value
        xshift = ground_altitude[Sim->XMAX-1][Sim->YMAX-1] - ground_altitude[0][Sim->YMAX-1];
        yshift = ground_altitude[Sim->XMAX-1][Sim->YMAX-1] - ground_altitude[Sim->XMAX-1][0];

        //y=starty;

        //#pragma omp parallel for default(none) private(x,y,x2,y2,nextx,nexty,xn,yn,lighttraveldist,xdist,ydist) shared(xsgn,ysgn,ax,ay,xyclaims,zslope,axdist,aydist)
        for (; x*xsgn < endx; x += xsgn) {
            for (y=starty; y*ysgn < endy; y += ysgn) {
                // when a new claim is made:
                // "try to shadow others"
                // Simple version, in which we assume that a tree has a vertically full shape (i.e. there is no gap inside in which light could pass) and there is no vertical overlap between two tress - so there can only be one single claim for each square.
                //step 1: check the next x'y' square in the opposite direction of the light source
                //step 2: evaluate a max_z that could be blocked the light (max_z corresponds to the first edge, closer to the sun)
                //step 3: if there is a claim located in this xy square, see if the altitude of the tree is above max_z
                //final steps: if so, then the x'y' claim of the corresponding tree is shadowed, do erase it (i.e. set v to -1); if not, then the current tree claim (x,y) as it received light

                if (xyclaims[x][y]->is_sunny != -10) {
                    //step 1: determine the corner that will likely cast the biggest shadow, in the current square
                    xn = 0; // the number of squares we are, along x axis
                    yn = 0; // .. .. .. .. .. .. .. .. .. .. .. .. y axis
                    x2 = ((double)x) + xsgn/2.0;
                    nextx = x2 + xsgn*ax/2.0;
                    y2 = ((double)y) + ysgn/2.0;
                    nexty = y2 + ysgn*ay/2.0;
                    lighttraveldist = 0;
                    //xdist = 0;
                    //ydist = 0;
                    current_altx = 0;
                    current_alty = 0;
                    //D printf("will start with %d,%d\n(reminder: ax=%f and ay=%f and abs(a)=%f)\n",x,y,ax,ay,fabs(a));
                    light_alt = 999999;
                    while (light_alt >= get_altitude(x+xn,y+yn)-1) { // stop condition is: does the lightbeam hit the ground at xn,yn ?
                        ////step 2: compute the coordinates of the first intersection with another square.
                        //D printf("at (%.2f,%.2f) [nextx=%.2f,nexty=%.2f]; next intersections are (%.2f,%.2f) and (%.2f,%.2f)\n",x2,y2,nextx,nexty,nextx,y+yn+ysgn,x+xn+xsgn,nexty);
                        if ( ysgn * (y+yn+ysgn - nexty) < 0 ) {
                            //D printf("intersected y=%f !\n",y+yn+ysgn);
                            lighttraveldist +=sqrt((x2-nextx)*(x2-nextx)+(y+yn+ysgn-y2)*(y+yn+ysgn-y2));
                            x2 = nextx;
                            nextx += ax*xsgn;
                            //nextx = x2+ax*xsgn;
                            y2 = y+yn+ysgn;
                            yn += ysgn;
                            //ydist += ysgn;
                        } else {
                            //D printf("intersected x=%f !\n",x+xn+xsgn);
                            lighttraveldist +=sqrt((y2-nexty)*(y2-nexty)+(x+xn+xsgn-x2)*(x+xn+xsgn-x2));
                            y2 = nexty;
                            nexty += ay*ysgn;
                            //nexty = y2+ay*ysgn;
                            x2 = x+xn+xsgn; //wtf
                            xn += xsgn;
                            //xdist += xsgn; 
                        }
                        //lighttraveldist = sqrt(xdist*xdist+ydist*ydist);
                        if (y+yn < 0) {
                            yn+=Sim->YMAX;
                            y2+=(double)Sim->YMAX;
                            nexty+=(double)Sim->YMAX;
                            current_alty--;
                        } else if (y+yn >= Sim->YMAX) {
                            yn-=Sim->YMAX;
                            y2-=(double)Sim->YMAX;
                            nexty-=(double)Sim->YMAX;
                            current_alty++;
                        }
                        if (x+xn < 0) {
                            xn+=Sim->XMAX;
                            x2+=(double)Sim->XMAX;
                            nextx+=(double)Sim->XMAX;
                            current_altx--;
                        } else if (x+xn >= Sim->XMAX) {
                            xn-=Sim->XMAX;
                            x2-=(double)Sim->XMAX;
                            nextx-=(double)Sim->XMAX;
                            current_altx++;
                        }
                        light_alt = xyclaims[x][y]->z - lighttraveldist*zslope - current_altx * xshift - current_alty * yshift;
                        // step 3: delete the claim if it is below the shadowed z
                        if (xyclaims[x+xn][y+yn]->z < light_alt) {
                            if ((x != x+xn) || (y != y+yn)) {
                                xyclaims[x+xn][y+yn]->is_sunny = 1;
                            } else {
                                printf("x=%d ; y=%d; xn=%d; yn=%d; dist=%f\n",x,y,xn,yn,lighttraveldist);
                            }
                            //break; // no need to continue to look if the ray hit an obstacle !
                        } // end if shadow claim
                    } // end while
                } // end if
            } // end y loop
        } // end x loop
    } else {
        // we reset xyclaims[x][y] to no interfer with further processing
        for (x = 0; x < Sim->XMAX; x++) {
            for (y = 0; y < Sim->YMAX; y++) {
                xyclaims[x][y]->is_sunny = 0;
            }
        }
    }
    //#pragma omp barrier

    // this loop over [x,y] updates the claims
    //#pragma omp parallel for default(none) private(x,y,t)
    for (x = 0; x < Sim->XMAX; x++) {
        for (y = 0; y < Sim->YMAX; y++) {
            if ( xyclaims[x][y]->z <= ground_altitude[x][y] || xyclaims[x][y]->is_sunny >= 1 ) {
                /* no trees around */
                biomass[y][x] = -1;
                profile[y][x] = 0;
            } else {
                //printf("we have a winner!\n");
                t = xyclaims[x][y]->t; /* the winner  */
                biomass[y][x] = t;
                profile[y][x] = (unsigned int) (xyclaims[x][y]->v * 100);
                tree[t]->ca++;
                if (Sim->SECTORS != 0) {
                    sector = xyclaims[x][y]->sec;
                    if ( tree[t]->cbhs[sector] > xyclaims[x][y]->v ) {
                        tree[t]->cbhs[sector] = xyclaims[x][y]->v;
                    }
                    if (tree[t]->ws[sector] < xyclaims[x][y]->rad) {
                        tree[t]->ws[sector] = xyclaims[x][y]->rad;
                    } else {
                        //printf("could not increase size to %f\n", xyclaims[x][y]->rad);
                    }
                } else {
                    tree[t]->cbh = xyclaims[x][y]->v;
                }
            } // end if trees around
        } // end y loop
    } // end x loop

    /*for (t = 0; t < stats[Sim->SMAX].tNow; t++) */
    //for (t = 0; t < num_trees; t++) {
    //	for (sector = 0; sector < Sim->SECTORS; sector++) {
    //		tree[t]->ws[sector] = sqrt(temptree[t]->ws[sector]); /* corrects for ws being in squared units */
    //		tree[t]->cbhs[sector] = temptree[t]->cbhs[sector]; 
    //	}
    //}


    /*now that you have realized ws, separate overstory and understory trees */
    for (t = 0; t < num_trees; t++) {
        if (tree[t]->ca <= 0) {
            tree[t]->is_ustory = 1; /* understory */
			tree[t]->mr = tree[t]->species->mrShade;    /* mortality rate for understory */
            for (sector = 0; sector < Sim->SECTORS; sector++) {
                tree[t]->sector_story[sector] = 1;
            }
        }
        else { /* (at least partly) overstory */
			tree[t]->mr = tree[t]->species->mrLight;    /* mortality rate for trees seeing light */
            for (sector = 0; sector < Sim->SECTORS; sector++) {
                if ( tree[t]->cbhs[sector] > 0.9 * tree[t]->h) {
                    // we put a cap on possible tree height at 90% of maximal height
                    tree[t]->cbhs[sector] = 0.9 * tree[t]->h;
                }
                if (tree[t]->ws[sector] < (tree[t]->max_sectors[sector] * 0.1)) {
                    /* this test aims to preserve a bit of an overstory crown
                     * which just got critically overtopped (loss > 90%) */
                    tree[t]->sector_story[sector] = 2; /* under and overstory */
                } else {
                    tree[t]->sector_story[sector] = 0; /*overstory */
                    /*since this has reached here, kill the understory */
                    tree[t]->ustory_ws[sector] = 0;
                    tree[t]->ustory_ca = 0;
                }
            }
            delta = 0;
            for (sector = 0; sector < Sim->SECTORS; sector++) {
                delta += tree[t]->sector_story[sector];
            }
            /*if the sum of the sectors is 0 then the tree is purely overstory */
            /*else it exists in both. */
            tree[t]->is_ustory = ((delta == 0) ? 0 : 2);
        }
//        if (tree[t]->is_ustory == 0) {
//            //stats[Sim->SMAX].tNow++;
//        } else {
//            //stats[Sim->SMAX].ustNow++;
//        }
    }
    // REWRITE CPP
    //std::cout << "FIXME: memory\n";
    //free(temptree);
}

/* UNDERSTORY GROWTH FUNCTION */
void grow_up_ustory() {
    int x, y, t;
    int num_trees = tree.size();
    short dx, dy; /* Used to determine periodic boundaries to find correct
                     sector. */
    /* dx, dy MUST be signed */
    // CPP
    //struct UNIT *temptree = (struct UNIT *)malloc(num_trees * sizeof(struct UNIT));
    //std::vector<struct UNIT *> temptree(num_trees);
    //for (int to_be_deleted_i=0; to_be_deleted_i<num_trees; to_be_deleted_i++) {
    //    temptree[to_be_deleted_i] =  (struct UNIT*) malloc(sizeof(struct UNIT));
    //}
    //temptree.reserve(num_trees);
    int sector;
    double sector_angle;
    double w2, r2, v, z;

    for (t = 0; t < num_trees; t++) {
        if ( tree[t]->is_ustory == 1 ) {
            // understory tree
            for (sector = 0; sector < Sim->SECTORS; sector++) {
                tree[t]->max_sectors[sector] = tree[t]->ustory_ws[sector];
                // the tree competes for understory -> let's dwarf its width and maximize its cbhs (they will be re-computed according to the xyuclaims)
                tree[t]->ws[sector] = 0.00;
                tree[t]->cbhs[sector] = tree[t]->h;
            }
        }
        //for (sector = 0; sector < Sim->SECTORS; sector++) {
        //	//temptree[t]->ws[sector] = 0.00;
        //	////temptree[t]->cbhs[sector] = tree[t]->cbhs[sector];
        //	//temptree[t]->cbhs[sector] = tree[t]->h;
        //}
        tree[t]->max_ws = max_ws(*tree[t]);
        tree[t]->ustory_ca = 0;
    }

    for (x = 0; x < Sim->XMAX; x++) {
        for (y = 0; y < Sim->YMAX; y++) {
            xyuclaims[x][y]->v = 0;
            xyuclaims[x][y]->z = ground_altitude[x][y];
            for (t = 0; t < num_trees; t++) {
                if ( tree[t]->is_ustory == 0) {
                    // overstory tree won't make under-story claims
                    continue;
                }

                /* distance */
                dx = x - tree[t]->x;
                dy = y - tree[t]->y;

                /* periodic boundaries: determines how the trees grow */
                if ((2 * abs(dx)) > Sim->XMAX) { /* if this goes over the edge of the picture */
                    if (dx > 0)
                        dx = dx - Sim->XMAX; /* dx positive */
                    else
                        dx = Sim->XMAX + dx;
                }
                if ((2 * abs(dy)) > Sim->YMAX) { /* same test as above */
                    if (dy > 0)
                        dy = dy - Sim->YMAX;
                    else
                        dy = Sim->YMAX + dy;
                }

                r2 = (double) ((dx * dx) + (dy * dy));

                w2 = tree[t]->max_ws * tree[t]->max_ws;
                if (r2 < w2) {   /* if it is a point within the projected crown */
                    if (Sim->SECTORS != 0) {
                        sector = get_sector(dx,dy);
                        sector_angle = (M_PI * 2.0 * ((double) sector) + 1.0) / ((double) Sim->SECTORS);
                        if (tree[t]->sector_story[sector] == 0) {
                            /*If this sector is overstory (only when the tree is both over- and under-story) */
                            continue;
                        }
                        w2 = 2.0+tree[t]->max_sectors[sector] * tree[t]->max_sectors[sector];
                    } else {
                        // handle the case where there are no sector
                        std::cerr << "no sector not implemented!\n" << std::endl;
                        sector = 0;
                        sector_angle = 0;
                    }
                    if (r2 < w2) {
                        v = tree[t]->ComputeClaim(r2, w2, sector);
                        z = v - sqrt(r2) * cos(sector_angle - Sim->CORIENTATION) * tan(Sim->CANGLE); /* the 2nd term corresponds to the change of altitude caused by the slope */
                        if (z<0 && v > 0) {
                            z = 0.001;
                        }
                        z += ground_altitude[x][y];

                        //if ( (RESTRICTCROWN == 0) || v >= tree[t]->cbhs[sector] ) {
                        if ( z > xyuclaims[x][y]->z || (z == xyuclaims[x][y]->z && (rand() % 2) ) ) {
                            xyuclaims[x][y]->rad = sqrt(r2);
                            xyuclaims[x][y]->t = t;
                            xyuclaims[x][y]->sec = sector;
                            xyuclaims[x][y]->v = v;
                            xyuclaims[x][y]->z = z;
                        }
                        //}
                    }
                }
            }
        }
    }

    //#pragma omp barrier

    //#pragma omp parallel for default(none) private(x,y,t)
    for (x = 0; x < Sim->XMAX; x++) {
        for (y = 0; y < Sim->YMAX; y++) {
            if ( xyuclaims[x][y]->z <= ground_altitude[x][y] ) {
                /* no trees around */
                ustory_biomass[y][x] = -1;
                ustory_profile[y][x] = 0;
            } else {
                t = xyuclaims[x][y]->t; /* the winner  */
                ustory_biomass[y][x] = t;
                ustory_profile[y][x] = (unsigned int) (xyuclaims[x][y]->v * 100);
                tree[t]->ustory_ca++;
                /*test to see if claim is greater than current cbh for sector */
                /*if so, replace. */
                if (Sim->SECTORS != 0) {
                    /* Below code checks EVERY POINT.
                     * If current point is farthest point in sector
                     * then you want to save it. The final longest radius
                     * is then found for the tree's sector's new radius.
                     */
                    sector = xyuclaims[x][y]->sec;
                    if ( tree[t]->cbhs[sector] > xyuclaims[x][y]->v ) {
                        tree[t]->cbhs[sector] = xyuclaims[x][y]->v;
                    }
                    if ( tree[t]->ws[sector] < xyuclaims[x][y]->rad) {
                        tree[t]->ws[sector] = xyuclaims[x][y]->rad;
                    }
                } else {
                    tree[t]->cbh = xyuclaims[x][y]->v;
                }
            }
        }
    }

    //for (t = 0; t < num_trees; t++) {
    //	if (tree[t]->is_ustory && tree[t]->ustory_ca >= 0) {
    //		for (sector = 0; sector < Sim->SECTORS; sector++) {
    //			//tree[t]->ustory_ws[sector] = sqrt(temptree[t]->ws[sector]); /* corrects for ws being in squared units */
    //			//if (tree[t]->cbhs[sector] < temptree[t]->cbhs[sector]) {
    //				// this additional check is intended for trees that are both under- and over-story
    //				//tree[t]->cbhs[sector] = temptree[t]->cbhs[sector]; 
    //			//}
    //		}
    //	}
    //}
}
/* END UNDERSTORY GROWTH FUNCTION */

void implicit_root_claim()
{
    int t;
    int num_trees = tree.size();
    double available_water = Sim->XMAX * Sim->YMAX * Sim->WATER_AVAILABLE;
    double requested_water = 0.0;
    for (t = 0; t < num_trees; t++) {
        tree[t]->root_w = 20.0 * tree[t]->d/2.0; // re-initialize root radius
        tree[t]->water_needed = (tree[t]->species->water_needed) * 3.14 * tree[t]->root_w * tree[t]->root_w;
        requested_water += tree[t]->water_needed;
    }
    double water_ratio = 1.0;
    std::cout << "available water=" <<available_water << "   requested water=" << requested_water << std::endl;
    if (requested_water > available_water) {
        water_ratio = available_water / requested_water;
    }
    for (t = 0; t < num_trees; t++) {
        tree[t]->water_uptake = tree[t]->water_needed * water_ratio;
    }
}

void circular_root_claim()
{
    int x, y, k, t;
    double w2, r2;
    int num_trees = tree.size();
    //short sector; // for future inclusion
    /* dx, dy MUST be signed */
    short dx, dy; /* Used to determine periodic boundaries to find correct
                     sector. */
    for (t = 0; t < num_trees; t++) {
        tree[t]->water_uptake = 0; // re-initialize water uptake
        tree[t]->root_w = 20.0 * tree[t]->d/2.0; // re-initialize root radius
        tree[t]->water_needed = (tree[t]->species->water_needed) * 3.14 * tree[t]->root_w * tree[t]->root_w;
        //tree[t]->water_needed = (0.5) * (1.0/3.0) * 3.14 * tree[t]->root_w * tree[t]->root_w;
        //tree[t]->water_needed = (tree[t]->species->water_needed) * (1.0/3.0) * 3.14 * tree[t]->root_w * tree[t]->root_w;
    }
    for (x = 0; x < Sim->XMAX; x++) {
        for (y = 0; y < Sim->YMAX; y++) {
            // 1> (re)initialize xyrc
            xyrc[x][y]->claim_nb = 0;
            xyrc[x][y]->norm_factor = 0;
            for (int k=0; k<CstMaxClaim; k++) {
                xyrc[x][y]->trees_i[k] = -1; // useless?
                xyrc[x][y]->trees_s[k] = 0; // useless?
            }
            // 2> add each tree to the list of xyrc, with the claim strength
            for (t = 0; t < num_trees; t++) {
                dx = x - tree[t]->x;
                dy = y - tree[t]->y;
                if ((2 * abs(dx)) > Sim->XMAX) {
                    if (dx > 0)
                        dx = dx - Sim->XMAX;
                    else
                        dx = Sim->XMAX + dx;
                }
                if ((2 * abs(dy)) > Sim->YMAX) {
                    if (dy > 0)
                        dy = dy - Sim->YMAX;
                    else
                        dy = Sim->YMAX + dy;
                }
                r2 = (double) ((dx * dx) + (dy * dy));
                w2 = tree[t]->root_w * tree[t]->root_w; /* squared root width */
                if (r2 < w2) { /* if point is within the root circle */
                    xyrc[x][y]->trees_i[xyrc[x][y]->claim_nb] = t;
                    xyrc[x][y]->trees_s[xyrc[x][y]->claim_nb] = 1.0 - sqrt(r2)/sqrt(w2);
                    //xyrc[x][y]->trees_s[xyrc[x][y]->claim_nb] = (tree[t]->d + 0.0001) / (1.0+sqrt(r2));
                    //xyrc[x][y]->trees_s[xyrc[x][y]->claim_nb] = 1.0/sqrt(r2);
                    xyrc[x][y]->norm_factor += xyrc[x][y]->trees_s[xyrc[x][y]->claim_nb]; // normalization factor
                    if (xyrc[x][y]->claim_nb < CstMaxClaim) {
                        xyrc[x][y]->claim_nb++;
                    } else {
                        std::cout << "Warning: one root claim was not taken into account!\n";
                    }
                    //}
                }
            }
        }
    }
    // 3> distribute the water
    for (x = 0; x < Sim->XMAX; x++) {
        for (y = 0; y < Sim->YMAX; y++) {
            for (k=0; k<xyrc[x][y]->claim_nb; k++) {
                if (xyrc[x][y]->norm_factor > 0) { // avoid NAN
                    tree[xyrc[x][y]->trees_i[k]]->water_uptake += water_content[x][y] * xyrc[x][y]->trees_s[k] / xyrc[x][y]->norm_factor;
                }
            }
        }
    }
}

//#ifdef __USE_ROOTS
//void water() {
//    int t;
//    const double water_demand = 100;/*, water_uptake = 10;*/
//    /*const double total_water = 10000;*/
//    int num_trees = stats[Sim->SMAX].tNow + stats[Sim->SMAX].ustNow;
//    /*long double tot_pot_up = 0;*/ /* total potential and realized uptakes */
//    /*double w_deficit = 0;*/ /* if total water < pot deficit = (pot - water)/pot */
//    /*double demand_coeff[3] = {1, 0.1, 0.5};*/ /*0:=overstory,1:=understory,2:=both */
//    /*double uptake_coeff[3] = {1, 2, 5};*/ /*0:=overstory,1:=understory,2:=both */
//    double ca_ratio = 0, water_ratio = 0;
//
//    for (t = 0; t < num_trees; t++) {
//        if (tree[t]->d >= 0 && tree[t]->d < 4)
//            tree[t]->wdemand = tree[t]->d * water_demand * species[tree[t]->s].waterdep[tree[t]->is_ustory];
//        else if (tree[t]->d >= 4)
//            tree[t]->wdemand = tree[t]->d * tree[t]->d * water_demand * species[tree[t]->s].waterdep[tree[t]->is_ustory];
//        /* tree[t]->wuptake = tree[t]->d * water_uptake * uptake_coeff[tree[t]->is_ustory];
//           tot_pot_up += tree[t]->wuptake; */
//    }
//    /* if (tot_pot_up > total_water)
//       w_deficit = (tot_pot_up - total_water)/tot_pot_up; */
//    for (t = 0; t < num_trees; t++) {
//        /* tree[t]->wrealuptake = (1 - w_deficit) * tree[t]->wuptake; */
//        ca_ratio = tree[t]->ca / tree[t]->pot_ca;
//        water_ratio = tree[t]->wrealuptake[ROOTMAX - 1] / tree[t]->wdemand;
//        /* printf("ca: %f, water: %f / %f\n", ca_ratio, tree[t]->wrealuptake[ROOTMAX - 1], tree[t]->wdemand); */
//        if (ca_ratio < water_ratio)
//            tree[t]->wlimit = 0;
//        else
//            tree[t]->wlimit = 1;
//    }	
//
//}

//int rootcmp(const void *a, const void *b) {
//    struct RCLAIM *a0, *b0;
//    a0 = (struct RCLAIM *)a;
//    b0 = (struct RCLAIM *)b;
//    if (a0->v > b0->v)
//        return 1;
//    if (a0->v < b0->v)
//        return -1;
//    return 0;
//}
//
//void root_competition(int *m, int *b) {
//    unsigned int t, sector, i, c = 0;
//    unsigned short x, y;
//    short dx, dy;
//    double w, r, max, sum;
//    /* Coefficients from http://www.wordiq.com/definition/Euclidean_distance */
//    const double e = 0.41, f = 0.941246;  
//    int num_trees = stats[Sim->SMAX].tNow + stats[Sim->SMAX].ustNow;
//    struct GROUND *tmproots = (struct GROUND *)malloc(num_trees * sizeof(struct GROUND));
//
//#if 0
//    for (t = 0; t < num_trees; t++) {
//        max = 0;
//        roots[t].total_top = 0;
//        for (sector = 0; sector < Sim->SECTORS; sector++) {
//            roots[t].rarea_top[sector] = 0;
//            tmproots[t].radii_top[sector] = 2.00;
//            if (roots[t].radii_top[sector] > max)
//                max = roots[t].radii_top[sector];
//        }
//        for (sector = 0; sector < ROOTMAX; sector++)
//            tree[t]->wrealuptake[sector] = 0;
//        tmproots[t].radius = max;
//        if (tree[t]->d >= midqualifier)
//            (*m)++;
//        else  
//            continue;
//        if (tree[t]->d >= botqualifier)
//            (*b)++;
//    }
//#endif
//
//    for (x = 0; x < Sim->XMAX; x++) {
//        for (y = 0; y < Sim->YMAX; y++) {
//            for (t = 0; t < num_trees; t++) {
//                if (tree[t]->d < 0) /* If tree is too small to compete at all */
//                    continue; /* Skip rest of iteration of loop */
//                if (!x && !y) {
//                    max = 0;
//                    roots[t].total_top = 0;
//                    for (sector = 0; sector < Sim->SECTORS; sector++) {
//                        roots[t].rarea_top[sector] = 0;
//                        tmproots[t].radii_top[sector] = 2.00;
//                        if (roots[t].radii_top[sector] > max)
//                            max = roots[t].radii_top[sector];
//                    }
//                    for (sector = 0; sector < ROOTMAX; sector++)
//                        tree[t]->wrealuptake[sector] = 0;
//                    tmproots[t].radius = max;
//                    if (tree[t]->d >= midqualifier)
//                        (*m)++;
//                    else  
//                        continue;
//                    if (tree[t]->d >= botqualifier)
//                        (*b)++;
//                }
//                /* distance */
//                dx = x - roots[t].x;
//                dy = y - roots[t].y;
//
//                /* toroidal bounds */
//                if ((2 * abs(dx)) > Sim->XMAX) {
//                    if (dx > 0)
//                        dx = dx - Sim->XMAX;
//                    else
//                        dx = Sim->XMAX + dx;
//                }
//                if ((2 * abs(dy)) > Sim->YMAX) {
//                    if (dy > 0)
//                        dy = dy - Sim->YMAX;
//                    else
//                        dy = Sim->YMAX + dy;
//                }
//
//                r = (double) ((dx * dx) + (dy * dy));
//
//                w = tmproots[t].radius * tmproots[t].radius;
//                if (r < w) { /* if it is a point within the super root area */
//                    if (Sim->SECTORS == 4) {
//                        if (dx > 0) {
//                            if (dy > 0)
//                                sector = 0;
//                            else
//                                sector = 3;
//                        }
//                        else {
//                            if (dy > 0)
//                                sector = 1;
//                            else
//                                sector = 2;
//                        }
//                    }	
//                    else if (Sim->SECTORS == 8) {	
//                        if (dx > 0) {
//                            if (dy > 0) {
//                                if (dx > dy)
//                                    sector = 0;
//                                else
//                                    sector = 1;
//                            }
//                            else {
//                                if (dx > -dy)
//                                    sector = 7;
//                                else
//                                    sector = 6;
//                            }
//                        }
//                        else {
//                            if (dy > 0) {
//                                if (-dx > dy)
//                                    sector = 3;
//                                else
//                                    sector = 2;
//                            }
//                            else {
//                                if (dx < dy)
//                                    sector = 4;
//                                else
//                                    sector = 5;
//                            }
//                        }
//                    }	
//                    /* sector is now the proper value */
//                    w = roots[t].radii_top[sector] * roots[t].radii_top[sector];
//                    if (r < w) { /* point is within that sector's area */	
//                        if (tree[t]->d >= 0) {
//                            rootclaims[c].t = t;
//                            /* calculate & save the (old) distance */
//                            /*	rootclaims[c].rad = exp((-0.5) * log(r)); */
//                            /* calculate the (new) distance & save it */
//                            if (dx > dy) {
//                                if (sector == 0)
//                                    rootclaims[c].rad = f * dx + e * dy;
//                                else if (sector == 7)
//                                    rootclaims[c].rad = f * dx - e * dy;
//                                else if (sector == 5)
//                                    rootclaims[c].rad = -1 * e * dx - f * dy;
//                                else if (sector == 6)
//                                    rootclaims[c].rad = e * dx - f * dy;
//                            }
//                            else {
//                                if (sector == 1)
//                                    rootclaims[c].rad = e * dx + f * dy;
//                                else if (sector == 2)
//                                    rootclaims[c].rad = f * dy - e * dx;
//                                else if (sector == 3)
//                                    rootclaims[c].rad = e * dy - f * dx;
//                                else if (sector == 4)
//                                    rootclaims[c].rad = -1 * e * dy - f * dx;
//                            }
//                            /* save the weight */
//                            rootclaims[c].v = tree[t]->d / rootclaims[c].rad;
//                            /* save the sector */
//                            rootclaims[c].sec = sector;
//                            /* save the actual radius^2 */
//                            rootclaims[c].rad2 = r;
//                            c++;
//                        }
//                    }
//                }	
//            }
//            if (c == 0) /* No roots around */
//                for (i = 0; i < ROOTMAX; i++) /* Spot is empty */
//                    roots_top[y][x][i] = Sim->TMAX;
//            else {
//                /*For those unfamiliar with it, qsort needs::
//                 * 1st: the array to sort
//                 * 2nd: # of elements in the array (cast to size_t)
//                 * 3rd: size of elems in array, e.g. an int arr => sizeof(int)
//                 * 4th: function pointer for comparison. MUST return an int
//                 *   pos: first element is greater than, 
//                 *   neg: first element is less than,
//                 *    0 : they're equal
//                 */
//                qsort(rootclaims, (size_t)c, sizeof(struct RCLAIM), rootcmp);
//                /* Let qsort sort the claims for me and use the last ROOTMAX
//                 * claims since c-1 is the largest, c-2 is the next largest,
//                 * etc. down to the last available slot.
//                 * -Ian Cordasco
//                 */
//                sum = 0;
//                /* we don't want to change c just yet, but we need the 
//                 * sum of the weights of the ROOTMAX winning trees */
//                for (i = 0, sector = c; (i < ROOTMAX) && (sector > 0); i++) 
//                    sum += rootclaims[--sector].v; /* We want to the sum of the largest claims */
//                for (i = 0; (i < ROOTMAX) && (c > 0); i++) {
//                    /* We want to fill the spot if we can, but if
//                     * c < ROOTMAX, we don't want junk values
//                     * so we want to make sure that c is also greater than zero
//                     */
//                    roots_top[y][x][i] = t = rootclaims[--c].t;
//                    sector = rootclaims[c].sec;
//                    roots[t].rarea_top[sector]++;
//                    roots[t].total_top++;
//                    if (tmproots[t].radii_top[sector] < rootclaims[c].rad2)
//                        tmproots[t].radii_top[sector] = rootclaims[c].rad2;
//                    w = top_water * rootclaims[c].v / sum;
//                    tree[t]->watersect[sector] += w;
//                    tree[t]->wrealuptake[0] += w;
//                    tree[t]->wrealuptake[ROOTMAX - 1] += w;
//                }
//                /* Now if c > ROOTMAX we have to make sure we start at the
//                 * beginning of the rootclaims array. If we don't we have 
//                 * the previous point's claims, c, minus ROOTMAX, and then the
//                 * next point's claims. This a) slows down the performance of
//                 * qsort since a significant portion of the array is already 
//                 * sorted, and b) is inaccurate since you're operating with the
//                 * previous point's claims. Now, that you're done with the claim
//                 * array, you can overwrite it from the beginning
//                 */
//                c = 0;
//            }
//        }
//    }
//    for (t = 0; t < num_trees; t++)
//        for (sector = 0; sector < Sim->SECTORS; sector++) {
//            roots[t].radii_top[sector] = sqrt(tmproots[t].radii_top[sector]);
//            /*Now that we know how well the sector has performed,
//             * find out it's ratio to the total and save that.*/
//            tree[t]->watersect[sector] /= tree[t]->wrealuptake[0];
//            if (isnan(tree[t]->watersect[sector])) 
//                /*If the ratio is so small that it is "NaN", set it equal to zero.*/
//                tree[t]->watersect[sector] = 0;
//        }
//
//    free(tmproots);
//}
//
//void rootmid_competition() {
//    unsigned int t, sector, i, c = 0;
//    unsigned short x, y;
//    short dx, dy;
//    double w, r, max, sum;
//    /* Coefficients from http://www.wordiq.com/definition/Euclidean_distance */
//    const double e = 0.41, f = 0.941246;  
//    int num_trees = stats[Sim->SMAX].tNow + stats[Sim->SMAX].ustNow;
//    struct GROUND *tmproots = (struct GROUND *)malloc(num_trees * sizeof(struct GROUND));
//
//#if 0
//    /* Reset or initialize some values crucial to the program */
//    for (t = 0; t < num_trees; t++) {
//        if (tree[t]->d >= midqualifier) {
//            max = 0;
//            roots[t].total_mid = 0;
//            for (sector = 0; sector < Sim->SECTORS; sector++) {
//                roots[t].rarea_mid[sector] = 0;
//                tmproots[t].radii_mid[sector] = 2.00;
//                /* Program is apparently sensitive to the temporary value */
//                if (roots[t].radii_mid[sector] > max)
//                    max = roots[t].radii_mid[sector];
//            }
//            tmproots[t].radius = max;
//        }
//    }
//#endif
//
//    /* Start competition for every point */
//    for (x = 0; x < Sim->XMAX; x++) {
//        for (y = 0; y < Sim->YMAX; y++) {
//            for (t = 0; t < num_trees; t++) {
//                if (tree[t]->d >= midqualifier) {
//                    /* NOTE: MOVE THIS OUTSIDE OF THIS if STATEMENT AFTER MERGE */
//                    if (!x && !y) {
//                        max = 0;
//                        roots[t].total_mid = 0;
//                        for (sector = 0; sector < Sim->SECTORS; sector++) {
//                            roots[t].rarea_mid[sector] = 0;
//                            tmproots[t].radii_mid[sector] = 2.00;
//                            /* Program is apparently sensitive to the temporary value */
//                            if (roots[t].radii_mid[sector] > max)
//                                max = roots[t].radii_mid[sector];
//                        }
//                        tmproots[t].radius = max;
//                    }
//                    /* If tree is too small to compete at all, it won't reach */
//                    /* distance */
//                    dx = x - roots[t].x;
//                    dy = y - roots[t].y;
//
//                    /* toroidal bounds */
//                    if ((2 * abs(dx)) > Sim->XMAX) {
//                        if (dx > 0)
//                            dx = dx - Sim->XMAX;
//                        else
//                            dx = Sim->XMAX + dx;
//                    }
//                    if ((2 * abs(dy)) > Sim->YMAX) {
//                        if (dy > 0)
//                            dy = dy - Sim->YMAX;
//                        else
//                            dy = Sim->YMAX + dy;
//                    }
//
//                    r = (double) ((dx * dx) + (dy * dy));
//
//                    w = tmproots[t].radius * tmproots[t].radius;
//                    if (r < w) { /* if it is a point within the super root area */
//                        if (Sim->SECTORS == 4) { 
//                            if (dx > 0) {
//                                if (dy > 0)
//                                    sector = 0;
//                                else
//                                    sector = 3;
//                            }
//                            else {
//                                if (dy > 0)
//                                    sector = 1;
//                                else
//                                    sector = 2;
//                            }
//                        }	
//                        else if (Sim->SECTORS == 8) {	
//                            if (dx >= 0) {
//                                if (dy > 0) {
//                                    if (dx > dy)
//                                        sector = 0;
//                                    else
//                                        sector = 1;
//                                }
//                                else {
//                                    if (dx > -dy)
//                                        sector = 7;
//                                    else
//                                        sector = 6;
//                                }
//                            }
//                            else {
//                                if (dy > 0) {
//                                    if (-dx > dy)
//                                        sector = 3;
//                                    else
//                                        sector = 2;
//                                }
//                                else {
//                                    if (dx < dy)
//                                        sector = 4;
//                                    else
//                                        sector = 5;
//                                }
//                            }
//                        }	
//                        /* sector is now the proper value */
//                        w = roots[t].radii_mid[sector] * roots[t].radii_mid[sector];
//                        if (r < w) { /* point is within that sector's area */	
//                            if (tree[t]->d >= midqualifier) {
//                                /* New competition uses weight: dbh/distance */
//                                /* save the (old) distance */
//                                /*	rootclaims[c].rad = exp((-0.5) * log(r)); */
//
//                                /* calculate the (new) distance & save it */
//                                rootclaims[c].t = t;
//                                /* Following assumes 8 sectors. 
//                                 * To add code for 4 sectors, look at link by
//                                 * the coefficients e and f at beginning of fn
//                                 * Also to better understand the logic go to
//                                 * that website.
//                                 */
//                                if (dx > dy) {
//                                    if (sector == 0)
//                                        rootclaims[c].rad = f * dx + e * dy;
//                                    else if (sector == 7)
//                                        rootclaims[c].rad = f * dx - e * dy;
//                                    else if (sector == 5)
//                                        rootclaims[c].rad = -1 * e * dx - f * dy;
//                                    else if (sector == 6)
//                                        rootclaims[c].rad = e * dx - f * dy;
//                                }
//                                else {
//                                    if (sector == 1)
//                                        rootclaims[c].rad = e * dx + f * dy;
//                                    else if (sector == 2)
//                                        rootclaims[c].rad = f * dy - e * dx;
//                                    else if (sector == 3)
//                                        rootclaims[c].rad = e * dy - f * dx;
//                                    else if (sector == 4)
//                                        rootclaims[c].rad = -1 * e * dy - f * dx;
//                                }
//                                /* save the weight for competition */
//                                rootclaims[c].v = tree[t]->d / rootclaims[c].rad;
//                                /* save the sector for later */
//                                rootclaims[c].sec = sector;
//                                /* save the actual radius^2 for realized set */
//                                rootclaims[c].rad2 = r;
//                                /* make sure you know you have claims */
//                                c++;
//                            }
//                        }
//                    }	
//                }
//            }
//            if (c == 0) /* No roots around */
//                for (i = 0; i < ROOTMAX; i++) /* Spot is empty */
//                    roots_mid[y][x][i] = Sim->TMAX;
//            else {
//                /*For those unfamiliar with it, qsort needs::
//                 * 1st: the array to sort
//                 * 2nd: # of elements in the array (cast to size_t)
//                 * 3rd: size of elems in array, e.g. an int arr => sizeof(int)
//                 * 4th: function pointer for comparison. MUST return an int
//                 *   pos: first element is greater than, 
//                 *   neg: first element is less than,
//                 *    0 : they're equal
//                 */
//                qsort(rootclaims, (size_t)c, sizeof(struct RCLAIM), rootcmp);
//                /* Let qsort sort the claims for me and use the last ROOTMAX
//                 * claims since c-1 is the largest, c-2 is the next largest,
//                 * etc. down to the last available slot.
//                 * -Ian Cordasco
//                 */
//                sum = 0;
//                /* we don't want to change c just yet, but we need the 
//                 * sum of the weights of the ROOTMAX winning trees */
//                for (i = 0, sector = c; (i < ROOTMAX) && (sector > 0); i++) 
//                    sum += rootclaims[--sector].v; /* We want to the sum of the largest claims */
//                for (i = 0; (i < ROOTMAX) && (c > 0); i++) {
//                    /* We want to fill the spot if we can, but if
//                     * c < ROOTMAX, we don't want junk values
//                     * so we want to make sure that c is also greater than zero
//                     */
//                    t = rootclaims[--c].t;
//                    roots_mid[y][x][i] = t;
//                    sector = rootclaims[c].sec;
//                    roots[t].total_mid++;
//                    if (tmproots[t].radii_mid[sector] < rootclaims[c].rad2)
//                        tmproots[t].radii_mid[sector] = rootclaims[c].rad2;
//                    w = mid_water * rootclaims[c].v / sum;
//                    if (rootclaims[c].rad2 == 0)
//                        for (dx = 0; dx < Sim->SECTORS; dx++)
//                            tree[t]->watersectmid[dx] += w/Sim->SECTORS;
//                    else{
//                        tree[t]->watersectmid[sector] += w;
//                        roots[t].rarea_mid[sector]++;
//                    }
//                    tree[t]->wrealuptake[1] += w;
//                    tree[t]->wrealuptake[ROOTMAX - 1] += w;
//                }
//                /* Now if c > ROOTMAX we have to make sure we start at the
//                 * beginning of the rootclaims array. If we don't we have 
//                 * the previous point's claims, c, minus ROOTMAX, and then the
//                 * next point's claims. This a) slows down the performance of
//                 * qsort since a significant portion of the array is already 
//                 * sorted, and b) is inaccurate since you're operating with the
//                 * previous point's claims. Now, that you're done with the claim
//                 * array, you can overwrite it from the beginning
//                 */
//                c = 0;
//            }
//        }
//    }
//    for (t = 0; t < num_trees; t++)
//        for (sector = 0; sector < Sim->SECTORS; sector++) {
//            if (tree[t]->d >= midqualifier) {
//                roots[t].radii_mid[sector] = sqrt(tmproots[t].radii_mid[sector]);
//                /*Now that we know how well the sector has performed,
//                 * find out it's ratio to the total and save that.*/
//                tree[t]->watersectmid[sector] /= tree[t]->wrealuptake[1];
//                if (isnan(tree[t]->watersectmid[sector])) {
//                    /*If the ratio is so small that it is "NaN", 
//                     *set it to an appropriate value.*/
//                    if (tree[t]->age < 15) {
//                        if (Sim->SECTORS == 8)
//                            /* If the tree is young give it 1/8th */
//                            tree[t]->watersectmid[sector] = 0.125;
//                        if (Sim->SECTORS == 4)
//                            /* If the tree is young give it 1/4th  */
//                            tree[t]->watersectmid[sector] = 0.25;
//                    }
//                    else
//                        /* Otherwise, we have no mercy */
//                        tree[t]->watersectmid[sector] = 0;
//                }
//            }
//        }
//
//    free(tmproots);
//}
//
//void rootbot_competition() {
//    unsigned int t, sector, i, c = 0;
//    unsigned short x, y;
//    short dx, dy;
//    double w, r, max, sum;
//    /* Coefficients from http://www.wordiq.com/definition/Euclidean_distance */
//    const double e = 0.41, f = 0.941246;  
//    int num_trees = stats[Sim->SMAX].tNow + stats[Sim->SMAX].ustNow;
//    struct GROUND *tmproots = (struct GROUND *)malloc(num_trees * sizeof(struct GROUND));
//
//#if 0
//    for (t = 0; t < num_trees; t++) {
//        if (tree[t]->d >= botqualifier) {
//            max = 0;
//            roots[t].total_mid = 0;
//            for (sector = 0; sector < Sim->SECTORS; sector++) {
//                roots[t].rarea_bot[sector] = 0;
//                tmproots[t].radii_bot[sector] = 2.00;
//                if (roots[t].radii_bot[sector] > max)
//                    max = roots[t].radii_bot[sector];
//            }
//            tmproots[t].radius = max;
//        }
//    }
//#endif
//
//    for (x = 0; x < Sim->XMAX; x++) {
//        for (y = 0; y < Sim->YMAX; y++) {
//            for (t = 0; t < num_trees; t++) {
//                if (tree[t]->d < botqualifier) /* If tree is too small to compete at all */
//                    continue; /* Skip rest of iteration of loop */
//                /* NOTE: MOVE ALL OF THESE TO BEFORE THE CONTINUE STATEMENTS */
//                if (!x && !y) {
//                    max = 0;
//                    roots[t].total_mid = 0;
//                    for (sector = 0; sector < Sim->SECTORS; sector++) {
//                        roots[t].rarea_bot[sector] = 0;
//                        tmproots[t].radii_bot[sector] = 2.00;
//                        if (roots[t].radii_bot[sector] > max)
//                            max = roots[t].radii_bot[sector];
//                    }
//                    tmproots[t].radius = max;
//                }
//                /* distance */
//                dx = x - roots[t].x;
//                dy = y - roots[t].y;
//
//                /* toroidal bounds */
//                if ((2 * abs(dx)) > Sim->XMAX) {
//                    if (dx > 0)
//                        dx = dx - Sim->XMAX;
//                    else
//                        dx = Sim->XMAX + dx;
//                }
//                if ((2 * abs(dy)) > Sim->YMAX) {
//                    if (dy > 0)
//                        dy = dy - Sim->YMAX;
//                    else
//                        dy = Sim->YMAX + dy;
//                }
//
//                r = (double) ((dx * dx) + (dy * dy));
//
//                w = tmproots[t].radius * tmproots[t].radius;
//                if (r < w) { /* if it is a point within the super root area */
//                    if (Sim->SECTORS == 4) {
//                        if (dx > 0) {
//                            if (dy > 0)
//                                sector = 0;
//                            else
//                                sector = 3;
//                        }
//                        else {
//                            if (dy > 0)
//                                sector = 1;
//                            else
//                                sector = 2;
//                        }
//                    }	
//                    else if (Sim->SECTORS == 8) {	
//                        if (dx > 0) {
//                            if (dy > 0) {
//                                if (dx > dy)
//                                    sector = 0;
//                                else
//                                    sector = 1;
//                            }
//                            else {
//                                if (dx > -dy)
//                                    sector = 7;
//                                else
//                                    sector = 6;
//                            }
//                        }
//                        else {
//                            if (dy > 0) {
//                                if (-dx > dy)
//                                    sector = 3;
//                                else
//                                    sector = 2;
//                            }
//                            else {
//                                if (dx < dy)
//                                    sector = 4;
//                                else
//                                    sector = 5;
//                            }
//                        }
//                    }	
//                    /* sector is now the proper value */
//                    w = roots[t].radii_bot[sector] * roots[t].radii_bot[sector];
//                    if (r < w) { /* point is within that sector's area */	
//                        if (tree[t]->d >= botqualifier) {
//                            /* New competition uses weight: dbh/distance */
//                            /* save the (old) distance */
//                            /*	rootclaims[c].rad = exp((-0.5) * log(r)); */
//                            /* calculate the (new) distance & save it */
//                            /* Following assumes 8 sectors. 
//                             * To add code for 4 sectors, look at link by
//                             * the coefficients e and f at beginning of fn
//                             * Also to better understand the logic go to
//                             * that website.
//                             */
//                            rootclaims[c].t = t;
//                            if (dx > dy) {
//                                if (sector == 0)
//                                    rootclaims[c].rad = f * dx + e * dy;
//                                else if (sector == 7)
//                                    rootclaims[c].rad = f * dx - e * dy;
//                                else if (sector == 5)
//                                    rootclaims[c].rad = -1 * e * dx - f * dy;
//                                else if (sector == 6)
//                                    rootclaims[c].rad = e * dx - f * dy;
//                            }
//                            else {
//                                if (sector == 1)
//                                    rootclaims[c].rad = e * dx + f * dy;
//                                else if (sector == 2)
//                                    rootclaims[c].rad = f * dy - e * dx;
//                                else if (sector == 3)
//                                    rootclaims[c].rad = e * dy - f * dx;
//                                else if (sector == 4)
//                                    rootclaims[c].rad = -1 * e * dy - f * dx;
//                            }
//                            /* save the weight */
//                            rootclaims[c].v = tree[t]->d / rootclaims[c].rad;
//                            /* save the sector */
//                            rootclaims[c].sec = sector;
//                            /* save the actual radius^2 */
//                            rootclaims[c].rad2 = r;
//                            c++;
//                        }
//                    }
//                }	
//            }
//            if (c == 0) /* No roots around */
//                for (i = 0; i < ROOTMAX; i++) /* Spot is empty */
//                    roots_bot[y][x][i] = Sim->TMAX;
//            else {
//                /*For those unfamiliar with it, qsort needs::
//                 * 1st: the array to sort
//                 * 2nd: # of elements in the array (cast to size_t)
//                 * 3rd: size of elems in array, e.g. an int arr => sizeof(int)
//                 * 4th: function pointer for comparison. MUST return an int
//                 *   pos: first element is greater than, 
//                 *   neg: first element is less than,
//                 *    0 : they're equal
//                 */
//                qsort(rootclaims, (size_t)c, sizeof(struct RCLAIM), rootcmp);
//                /* Let qsort sort the claims for me and use the last ROOTMAX
//                 * claims since c-1 is the largest, c-2 is the next largest,
//                 * etc. down to the last available slot.
//                 * -Ian Cordasco
//                 */
//                sum = 0;
//                /* we don't want to change c just yet, but we need the 
//                 * sum of the weights of the ROOTMAX winning trees */
//                for (i = 0, sector = c; (i < ROOTMAX) && (sector > 0); i++) 
//                    sum += rootclaims[--sector].v; /* We want to the sum of the largest claims */
//                for (i = 0; (i < ROOTMAX) && (c > 0); i++) {
//                    /* We want to fill the spot if we can, but if
//                     * c < ROOTMAX, we don't want junk values
//                     * so we want to make sure that c is also greater than zero
//                     */
//                    t = rootclaims[--c].t;
//                    roots_bot[y][x][i] = t;
//                    sector = rootclaims[c].sec;
//                    roots[t].total_bot++;
//                    if (tmproots[t].radii_bot[sector] < rootclaims[c].rad2)
//                        tmproots[t].radii_bot[sector] = rootclaims[c].rad2;
//                    w = bot_water * rootclaims[c].v / sum;
//                    if (rootclaims[c].rad2 == 0)
//                        for (dx = 0; dx < Sim->SECTORS; dx++)
//                            tree[t]->watersectbot[dx] += w/Sim->SECTORS;
//                    else{
//                        tree[t]->watersectbot[sector] += w;
//                        roots[t].rarea_bot[sector]++;
//                    }
//                    tree[t]->wrealuptake[2] += w;
//                    tree[t]->wrealuptake[ROOTMAX - 1] += w;
//                }
//                /* Now if c > ROOTMAX we have to make sure we start at the
//                 * beginning of the rootclaims array. If we don't we have 
//                 * the previous point's claims, c, minus ROOTMAX, and then the
//                 * next point's claims. This a) slows down the performance of
//                 * qsort since a significant portion of the array is already 
//                 * sorted, and b) is inaccurate since you're operating with the
//                 * previous point's claims. Now, that you're done with the claim
//                 * array, you can overwrite it from the beginning
//                 */
//                c = 0;
//            }
//        }
//    }
//    for (t = 0; t < num_trees; t++)
//        for (sector = 0; sector < Sim->SECTORS; sector++) {
//            if (tree[t]->d >= botqualifier) {
//                /* take the realized radius for each sector */
//                roots[t].radii_bot[sector] = sqrt(tmproots[t].radii_bot[sector]);
//                /*Now that we know how well the sector has performed,
//                 * find out it's ratio to the total and save that.*/
//                tree[t]->watersectbot[sector] /= tree[t]->wrealuptake[1];
//                if (isnan(tree[t]->watersectbot[sector])) {
//                    /*If the ratio is so small that it is "NaN", 
//                     * set it to something.*/
//                    if (tree[t]->age < 100) {
//                        if (Sim->SECTORS == 8)
//                            tree[t]->watersectbot[sector] = 0.125;
//                        if (Sim->SECTORS == 4)
//                            tree[t]->watersectbot[sector] = 0.25;
//                    }
//                    else
//                        /* Otherwise we have no mercy */
//                        tree[t]->watersectbot[sector] = 0;
//                }
//            }
//        }
//    free(tmproots);
//}
//#endif

void sunlight() {
    double dx, dy, sumx_neg, sumx_pos, sumy_neg, sumy_pos; /*new movement */
    double growth_ratio, growth_ratio2; 
#ifdef __USE_ROOTS
//    double water_ratio = 0;
#endif
    int n;
    unsigned int t;
    unsigned int num_trees = tree.size();
    double increment_factor;

    int todeletecountofrestricted = 0;

    for (t = 0; t < num_trees; t++) {
        // another way: (untested)
        //growth_ratio = tree[t]->tr + (1 - tree[t]->tr) * tree[t]->ca / tree[t]->pot_ca;
        // usual way:
        growth_ratio = tree[t]->ca / tree[t]->pot_ca;
        if (growth_ratio < tree[t]->tr) {
            growth_ratio = tree[t]->tr;
        }
        if (Sim->ROOT >= 1) {
            if (tree[t]->water_uptake < tree[t]->water_needed) {
                growth_ratio2 = tree[t]->water_uptake / tree[t]->water_needed;
                if (Sim->LIGHTWATERGROWTH == 0) { // use only sunlight
                    growth_ratio = growth_ratio; // silly but explicit
                } else if (Sim->LIGHTWATERGROWTH == 1) { // multiplicative influence
                    growth_ratio = growth_ratio * growth_ratio2;
                } else { // LIGHTWATERGROWTH==2 => law of the minimum
                    if (growth_ratio2 < growth_ratio) {
                        growth_ratio = growth_ratio2;
                    }
                }
                todeletecountofrestricted++;
            }
        }
        if (tree[t]->is_ustory == 0) { 
            /*tree[t]->d += diameter_inc(tree[t]); */
            // // // OLD VERSION tree[t]->d += diameter_inc(tree[t]) * ca_ratio; 
            // // // not suitable, because it leads to a decrease of trunk diameter growth as soon as a tree reaches canopy level
            tree[t]->d += diameter_inc(*tree[t]) * growth_ratio; // re-enabled after correction for minimal ca_ratio=tree[t]->tr
            // tree[t]->d += diameter_inc(*tree[t]) * (ca_ratio + tree[t]->tr); // CPP
            /*if ((tree[t]->ca / pre_crown(tree[t])) > tree[t]->tr)    trees having some CONSIDERABLE canopy */
            /*tree[t]->d += diameter_inc(tree[t]) * tree[t]->ca / pre_crown(tree[t]);  */
        } else { /* trees ALMOST completely in shadows -- understory */
            tree[t]->d += diameter_inc(*tree[t]) * tree[t]->tr;
        }
        if (tree[t]->d > tree[t]->species->dMax) { // May 31, enforce maximal diameter to prevent deviation caused by PRI
            //std::cout << "tree with diameter: " << tree[t]->d << " but max diameter is " << tree[t]->species->dMax << std::endl;
            //std::cout << "otherwise, mr is " << tree[t]->species->mrShade << std::endl;
            tree[t]->d = tree[t]->species->dMax;
        }
        /* new tree sizes */
        if (Sim->VERTICALGROWTH == 0) {
            /* with an increment proportional to the trunk diameter: */
            //printf("tree %d, old h=%.2f, new h=%.2f\n",t,tree[t]->h,height(tree[t]));
            tree[t]->h = height(*tree[t]);
        } else if (Sim->VERTICALGROWTH == 1) {
            /* with an increment proportional to the vertical sunlit canopy: */
            increment_factor = 0;
            for (n = 0; n < Sim->SECTORS; n++) {
                // store the mean difference between the height and the cbh of the sectors
                increment_factor += (tree[t]->h - tree[t]->cbhs[n])/Sim->SECTORS;
            }
            increment_factor = increment_factor / tree[t]->h;
            if (increment_factor < 0.5) {
                increment_factor = 0.5;
            }
            if (Sim->RESTRICTHEIGHT == 0) {
                /* without limitation: */
                //increment_factor *= 1.0;
            } else if (Sim->RESTRICTHEIGHT == 1) {
                /* with a cap on maximal radius: */
                if (tree[t]->h < tree[t]->species->hMax) {
                    //increment_factor *= 1.0;
                } else {
                    increment_factor = 0;
                }
            } else if (Sim->RESTRICTHEIGHT == 2) {
                /* with an increment proportional to the ratio of current/max heights: */
                increment_factor *= (1 - tree[t]->h / tree[t]->species->hMax);
            }

            if (increment_factor < 0) {
                /* to prevent a particularly nasty bug in which the crown gets inversed  */
                increment_factor = 0;
            } else if (increment_factor > 1) {
                increment_factor = 1;
            }
            tree[t]->h += tree[t]->species->hinc * increment_factor;
        }
        /* tree[t]->w = pre_radius(tree[t]); */

        /* calculate each sector's new radius  */
        for (n = 0; n < Sim->SECTORS; n++) {
            if (tree[t]->sector_story[n] == 0) { 
                /* canopy */
                if (Sim->RESTRICTRADIUS == 0) {
                    /* without limitation: */
                    increment_factor = 1.0;
                } else if (Sim->RESTRICTRADIUS == 1) {
                    /* with a cap on maximal radius: */
                    if (tree[t]->ws[n] < tree[t]->species->wMax) {
                        increment_factor = 1.0;
                    } else {
                        increment_factor = 0.0;
                    }
                } else if (Sim->RESTRICTRADIUS == 2) {
                    /* with an increment proportional to the ratio of current/max widths: */
                    increment_factor = (1.0 - tree[t]->ws[n] / tree[t]->species->wMax);
                } else if (Sim->RESTRICTRADIUS == 3) {
                    /* with sophisticated limitation:
                     * we take here the inverse formula of the shape (used previously to compute the vertical claim):
                     * new_radius = species_increment * (1 - old_radius / maximal_radius_at_given_height)
                     * with:
                     * maximal_radius_at_given_height = maximal_radius * (1 - current_height/maximal_height)^species_exp_constant */
                    increment_factor = 1.0*(1.0 - (tree[t]->ws[n] / (tree[t]->species->wMax * exp(log(1.0 - (tree[t]->h/tree[t]->species->hMax))*tree[t]->species->csexp)) ));
                } else {
                    std::cerr << "\nInvalid value for RESTRICTRADIUS\nEXITING!!!\n" << std::endl;
                    exit(1);
                }
                if (increment_factor < 0.0) {
                    /* to prevent a particularly nasty bug in which the crown gets inversed  */
                    increment_factor = 0.0;
                } else if (increment_factor > 1.0) {
                    increment_factor = 1.0;
                }
                // tree[t]->ws[n] += tree[t]->species->rinc * increment_factor; // old rinc
                tree[t]->ws[n] += crown_inc(*tree[t]) * increment_factor;
                if (Sim->POTENTIALCROWN == 2) {
                    tree[t]->pot_radius += crown_inc(*tree[t]) * increment_factor;
                }
            } else {
                /* understory */
                if (Sim->RESTRICTRADIUS == 0) {
                    /* without limitation: */
                    increment_factor = 1;
                } else if (Sim->RESTRICTRADIUS == 1) {
                    /* with a cap on maximal radius: */
                    if (tree[t]->ustory_ws[n] < tree[t]->species->wMax) {
                        increment_factor = 1;
                    } else {
                        increment_factor = 0;
                    }
                } else if (Sim->RESTRICTRADIUS == 2) {
                    /* with an increment proportional to the ratio of current/max widths: */
                    increment_factor = (1 - tree[t]->ustory_ws[n] / tree[t]->species->wMax);
                } else if (Sim->RESTRICTRADIUS == 3) {
                    /* with sophisticated limitation:
                     * we take here the inverse formula of the shape (used previously to compute the vertical claim):
                     * new_radius = species_increment * (1 - old_radius / maximal_radius_at_given_height)
                     * with:
                     * maximal_radius_at_given_height = maximal_radius * (1 - current_height/maximal_height)^species_exp_constant */
                    increment_factor = 1*(1 - (tree[t]->ustory_ws[n] / (tree[t]->species->wMax * exp(log(1 - (tree[t]->h/tree[t]->species->hMax))*tree[t]->species->csexp)) ));
                } else {
                    std::cerr << "\nInvalid value for RESTRICTRADIUS\nEXITING!!!\n" << std::endl;
                    exit(1);
                }
                if (increment_factor < 0) {
                    /* to prevent a particularly nasty bug in which the crown gets inversed  */
                    increment_factor = 0;
                } else if (increment_factor > 1) {
                    increment_factor = 1;
                }
                if (tree[t]->sector_story[n] == 1) {
                    /* understory */
                    tree[t]->ustory_ws[n] += tree[t]->species->rinc * tree[t]->tr * 3.0 * increment_factor;
                    // original version: the factor is 3
                } else {
                    /* under and over-story */
                    tree[t]->ustory_ws[n] += tree[t]->species->rinc * tree[t]->tr * 4.0 * increment_factor;
                    // original version: the factor is 4
                }
            }
        }
        /* MOVEMENT */
        if (Sim->SECTORS == 8) {
            /* Along the abscissa axis */
            sumx_neg = tree[t]->ws[2] + tree[t]->ws[3] + tree[t]->ws[4] + tree[t]->ws[5];
            /*sumx_neg is the left half of the abscissa */
            sumx_pos = tree[t]->ws[0] + tree[t]->ws[1] + tree[t]->ws[6] + tree[t]->ws[7];
            /*sumx_pos is the right half of the abscissa */
            /* Along the ordinate axis */
            sumy_neg = tree[t]->ws[4] + tree[t]->ws[5] + tree[t]->ws[6] + tree[t]->ws[7];
            /*sumy_neg is the bottom half of the ordinate */
            sumy_pos = tree[t]->ws[0] + tree[t]->ws[1] + tree[t]->ws[2] + tree[t]->ws[3];
            /*sumy_pos is the top half of the ordinate */
        }
        else if (Sim->SECTORS == 4) {
            /* Along the abscissa axis */
            /*sumx_neg is the left half of the abscissa */
            sumx_neg = tree[t]->ws[1] + tree[t]->ws[2];
            /*sumx_pos is the right half of the abscissa */
            sumx_pos = tree[t]->ws[0] + tree[t]->ws[3];
            /* Along the ordinate axis */
            /*sumy_neg is the bottom half of the ordinate */
            sumy_neg = tree[t]->ws[2] + tree[t]->ws[3];
            /*sumy_pos is the top half of the ordinate */
            sumy_pos = tree[t]->ws[0] + tree[t]->ws[1];
        } else {
            sumx_neg = 0;
            sumx_pos = 0;
            sumy_neg = 0;
            sumy_pos = 0;
        }
        if (sumx_pos - sumx_neg > 0) {
            dx = (tree[t]->h * (1 - (sumx_neg/sumx_pos)) * tree[t]->species->tiltMax);
            if (dx < 1)
                dx = 0;
            else if (dx > 2)
                dx = 1;
            dx += (double)tree[t]->x;
            if (dx > (double)Sim->XMAX)
                dx = dx - (double)Sim->XMAX;
            tree[t]->x = (short)dx;
        }
        else if (sumx_pos - sumx_neg < 0) {
            dx = (tree[t]->h * (1 - (sumx_pos/sumx_neg)) * tree[t]->species->tiltMax);
            if (dx < 1)
                dx = 0;
            else if (dx > 2)
                dx = 1;
            dx = (double)tree[t]->x - dx;
            if (dx < 0)
                dx = Sim->XMAX + dx;
            tree[t]->x = (short)dx;
        }
        if (sumy_pos - sumy_neg > 0) {
            dy = (tree[t]->h * (1 - (sumy_neg/sumy_pos)) * tree[t]->species->tiltMax);
            if (dy < 1)
                dy = 0;
            else if (dy > 2)
                dy = 1;
            dy += (double)tree[t]->y;
            if (dy > (double)Sim->YMAX)
                dy = dy - (double)Sim->YMAX;
            tree[t]->y = (short)dy;
        }
        else if (sumy_pos - sumy_neg < 0) {
            dy = (tree[t]->h * (1 - (sumy_pos/sumy_neg)) * tree[t]->species->tiltMax);
            if (dy < 1)
                dy = 0;
            else if (dy > 2)
                dy = 1;
            dy = (double)tree[t]->y - dy;
            if (dy < 0)
                dy = Sim->YMAX + dy;
            tree[t]->y = (short)dy;
        }
    }
    std::cout << "\n    " << todeletecountofrestricted << " trees limited over a total pop of " << tree.size() << "\n";
    //std::cout << "\n   tree 10 water uptake: " << tree[10]->water_uptake << " VS tree 10 water needed: " << tree[10]->water_needed << "\n";
}

//void measure() { /* crown shape analysis */
//	int s, t, n;
//	unsigned short x, y;
//	int num_trees = tree.size();
//	short xx, yy;
//	/*double h, hMax, hMin, area; */
//	double prof_y_x;
//
//	for (s = 0; s <= Sim->SMAX; s++) {    /* reset species statistics data */
//		stats[s].ca = 0;
//		stats[s].baCan = 0;
//		stats[s].baAll = 0;
//		stats[s].tActive = 0;
//	}
//
//	for (t = 0; t < num_trees; t++) {
//		stats[tree[t]->s].baAll += tree[t]->d * tree[t]->d * M_PI / 4;    /* sum of basal area of all trees of one species */
//		stats[Sim->SMAX].baAll += tree[t]->d * tree[t]->d * M_PI / 4; /* sum of dbh of all trees of all species */
//
//		//if (tree[t]->ca > 0) {
//		if (tree[t]->is_ustory == 0) { // If overstory
//			stats[tree[t]->s].tActive++;
//			stats[Sim->SMAX].tActive++;
//			stats[tree[t]->s].ca += tree[t]->ca;   /* sum of canopy area of one cpecies */
//			stats[Sim->SMAX].ca += tree[t]->ca;    /* sum of canopy area of all cpecies */
//			stats[tree[t]->s].baCan += tree[t]->d * tree[t]->d * M_PI / 4;    /* sum of dbh of canopy trees of one species */
//			stats[Sim->SMAX].baCan += tree[t]->d * tree[t]->d * M_PI / 4; /* sum of dbh of canopy trees of all species */
//			tree[t]->mr = species[tree[t]->s].mrLight;    /* mortality rate for trees seeing light */
//		}
//		else
//			tree[t]->mr = species[tree[t]->s].mrShade;    /* mortality rate for trees in the shadow */
//		/*if (tree[t]->wlimit == 1) */
//		/*    tree[t]->mr *= (tree[t]->wdemand/tree[t]->wrealuptake); */
//		tree[t]->cbh = tree[t]->h;
//		tree[t]->brdh = 0;
//		tree[t]->brdl = 0;
//		tree[t]->water_mr = 0.01; // fixed additional mortality rate of water-limited trees
//#ifdef WATER_GRADIENT
//		tree[t]->water_mr = 0.2; // huge mortality to show transition to grassland
//#endif
//#ifdef ESA_WATER
//		tree[t]->water_mr = 0.2; // huge mortality to show transition to grassland
//#endif
//#ifdef LOWER_DIST
//        tree[t]->mr *= 0.5;
//#endif
//#ifdef HIGHER_DIST
//        tree[t]->mr *= 2;
//#endif
//
//	}
//
//	for (s = 0; s < Sim->SMAX; s++) { /* set number of newborns for the next cycle */
//		/*      species[s].nSeed=(short)((double)Sim->nSeedTotal*stats[s].baAll/stats[Sim->SMAX].baAll); */ /*number of newborns ~ sum of ba */ 
//		/*      species[s].nSeed=(short)((double)Sim->nSeedTotal*stats[s].baCan/stats[Sim->SMAX].baCan);*/    /* number of newborns ~ sum of canopy ba */
//
//		species[s].nSeed = (short) ((double) Sim->nSeedTotal * stats[s].ca / stats[Sim->SMAX].ca);  /* number of newborns ~ sum of canopy */
//	}
//
//	/* crown border length and average height */
//
//	for (x = 0; x < Sim->XMAX; x++) {
//		for (y = 0; y < Sim->YMAX; y++) {
//			n = 0;
//			t = biomass[y][x];
//			prof_y_x = ((double) profile[y][x] / 100);
//			if (t != -1) {
//				/* set the crown base height as the lowest point of the current crown */
//				if (tree[t]->cbh > prof_y_x) {
//					tree[t]->cbh = prof_y_x;
//                }
//				/* calculate number of neighbour leaves */
//				for (yy = -1; yy <= 1; yy++) {
//					for (xx = -1; xx <= 1; xx++) {
//						if (biomass[(unsigned) (y + yy) % Sim->YMAX][(unsigned) (x + xx) % Sim->XMAX] == t) {
//							n++;
//                        }
//                    }
//                }
//			}
//			if (n > 3 && n < 8) {    /* crown border criterium */
//				tree[t]->brdl++; /* number of points ~ crown border length */
//				tree[t]->brdh += prof_y_x;   /* crown border point height in dm */
//			}
//		} 
//    }
//    for (t = 0; t < num_trees; t++) {
//		if (tree[t]->brdl > 0) {
//            /* tree crown has a border */
//			tree[t]->brdh /= (double) tree[t]->brdl;  /* average crown border height in dm */
//        }
//    }
//
//	printf("Total: %4i(%3i)\n", num_trees, stats[Sim->SMAX].tActive);
//	}

inline void set_dxdy(short int * dx, short int * dy) {
    /* periodic boundaries: determines how the trees grow */
    if ((2 * abs(*dx)) > Sim->XMAX) { /* if this goes over the edge of the picture */
        if (*dx > 0)
            *dx = *dx - Sim->XMAX; /* dx positive */
        /* If the point goes off the positive side of
         * the plot, then it becomes negative since it
         * is then placed to the left of the tree's
         * root.
         * Else, it is negative and becomes positive.
         * The same logic applies to dy.
         */
        else
            *dx = Sim->XMAX + *dx;
        /* Since dx is negative, adding it to Sim->XMAX is
         * equivalent to Sim->XMAX - abs(dx). This is just
         * faster since you're not first calling 
         * another function.
         * Same applies to dy.
         */
    }
    if ((2 * abs(*dy)) > Sim->YMAX) { /* same test as above */
        if (dy > 0)
            *dy = *dy - Sim->YMAX;
        else
            *dy = Sim->YMAX + *dy;
    }
}

inline int get_sector(int dx, int dy) {
    int sector;
    /* SECTOR LOGIC:
     * For 4 sectors:
     *   If dx & dy are positive, then in quadrant 1,
     *   If dx is negative & dy is positive, then in quad 2,
     *   If dx & dy are negative, then in quadrant 3,
     *   If dx is positive & dy is negative, then in quad 4.
     * For 8 sectors:
     *   If dx & dy are positive:
     *     If dx > dy: 1st octant [0, pi/4]
     *     If dy > dx: 2nd octant [pi/4, pi/2]
     *   If dx < 0, dy > 0:
     *     If dy > |dx|: 3rd octant [pi/2, 3pi/4]
     *     If |dx| > dy: 4th octant [3pi/4, pi]
     *   If dx & dy < 0:
     *     If |dx| > |dy|: 5th octant [pi, 5pi/4]
     *     If |dy| > |dx|: 6th octant [5pi/4, 3pi/2]
     *   If dx > 0, dy < 0:
     *     If |dy| > dx: 7th octant [3pi/2, 7pi/2]
     *     If dx > |dy|: 8th octant [7pi/2, 2pi]
     * -Ian Cordasco <icordasc@math.stevens.edu> Summer 2010
     */
    if (Sim->SECTORS == 4) {
        if (dx > 0) {
            if (dy > 0) { /* Quadrant: 1 */
                sector = 0;
            }
            else { /* Quadrant: 4 */
                sector = 3;
            }
        }
        else {
            if (dy > 0) { /* Quadrant: 2 */
                sector = 1;
            }
            else { /* Quadrant: 3 */
                sector = 2;
            }
        }
    }
    else if (Sim->SECTORS == 8) {
        if (dx > 0) { 
            if (dy > 0) {
                if (dx > dy) { /* Octant: 1 */
                    sector = 0;
                }
                else { /* Octant: 2 */
                    sector = 1;
                }
            }
            else { 
                if (dx > -dy) { /* Octant: 8 */
                    sector = 7;
                }
                else { /* Octant: 7 */
                    sector = 6;
                }
            }
        }
        else {
            if (dy > 0) { 
                if (-dx > dy) { /* Octant: 4 */
                    sector = 3;
                }
                else { /* Octant: 3 */
                    sector = 2;
                }
            }
            else { 
                if (dx < dy) { /*If dx is more negative than dy */
                    /* Octant: 5 */
                    sector = 4;
                }
                else { 
                    /* Octant: 6 */
                    sector = 5;
                }
            }
        }
    } else {
        sector = 0;
    }
    return sector;
}



void color_bitmap_mainroot(void) {       
    /* color the crown */
    int x, y, k;
    int winner; // the main sucker
    double water_uptake_xy; // max water uptake at one point of space
    for (x = 0; x < Sim->XMAX; x++) {
        for (y = 0; y < Sim->YMAX; y++) {
            winner = -1;
            water_uptake_xy = 0;
            for (k=0; k<xyrc[x][y]->claim_nb; k++) {
                if (xyrc[x][y]->trees_s[k] > water_uptake_xy) {
                    winner = xyrc[x][y]->trees_i[k];
                    water_uptake_xy = xyrc[x][y]->trees_s[k];
                }
            }
            if (winner == -1) {
                bitmap[y][x] = pBlack;
            } else {
                // bitmap[y][x] = tree[winner]->color;
                bitmap[y][x].B =
                    (unsigned char) (tree[winner]->color.B * (water_uptake_xy / xyrc[x][y]->norm_factor));
                bitmap[y][x].G =
                    (unsigned char) (tree[winner]->color.G * (water_uptake_xy / xyrc[x][y]->norm_factor));
                bitmap[y][x].R =
                    (unsigned char) (tree[winner]->color.R * (water_uptake_xy / xyrc[x][y]->norm_factor));
                //bitmap[y][x] = pGreen;
            }
        }
    }
    //std::cout << "gradient at 252x252: " << xyrc[252][252]->trees_s[winner] / xyrc[252][252]->norm_factor << "\n";
}

void color_bitmap_cnp(void) {       
    /* color the crown */
    int x, y, t;
    int num_trees = tree.size();
    for (x = 0; x < Sim->XMAX; x++)
        for (y = 0; y < Sim->YMAX; y++) {
            /* canopy distribution */
            if (biomass[y][x] == -1) {
                bitmap[y][x] = pBlack;
            } else {
                bitmap[y][x] = tree[biomass[y][x]]->color;
            }
        }
    /* black dot at the center of growth */
    for (t = 0; t < num_trees; t++)
        if (!(tree[t]->is_ustory)) {
            bitmap[tree[t]->y][tree[t]->x].R = 0;
            bitmap[tree[t]->y][tree[t]->x].G = 0;
            bitmap[tree[t]->y][tree[t]->x].B = 0;
        }
}

void color_bitmap_hgt(void)
{
    int x, y, t;
    int num_trees = tree.size();

    double gradient;

    /* color the crown */
    for (x = 0; x < Sim->XMAX; x++) {
        for (y = 0; y < Sim->YMAX; y++) {
            if (biomass[y][x] == -1) {
                bitmap[y][x].R = 0x00;
                bitmap[y][x].G = 0x00;
                bitmap[y][x].B = 0x00;
            } else {
                /* height distribution */
                gradient = ((double)profile[y][x]) / (100.0 * ((double)Sim->HSCALE) * (1.0 + hRnd));
                bitmap[y][x].B = (unsigned char) (tree[biomass[y][x]]->species->leaf.B * gradient);
                bitmap[y][x].G = (unsigned char) (tree[biomass[y][x]]->species->leaf.G * gradient);
                bitmap[y][x].R = (unsigned char) (tree[biomass[y][x]]->species->leaf.R * gradient);
            }

                // // removed this part:
                // //                 /* crown border points */
                // //                 int n=0;
                // //                 for (int yy = -1; yy <= 1; yy++)
                // //                     for (int xx = -1; xx <= 1; xx++)
                // //                         if (biomass[(y + yy) % Sim->YMAX][(x + xx) % Sim->XMAX] == biomass[y][x])
                // //                             n++;
                // //             }
                // //             if (n > 3 && n < 8)  /* crown border criterium */
                // //                 if (Sim->SHOWTREELIMITS)
                // //                     bitmap[y][x] = pYellow; /* yellow line at the crown border */
        }
    }


    // // Code to show the roots as well
    for (t = 0; t < num_trees; t++) {
        /*if (tree[t]->ca > 0) { */   /* tree has a canopy */
        if ((tree[t]->is_ustory == 0) || (tree[t]->is_ustory == 2)) {
            if (Sim->SHOWTREEROOT) {  /* black cross at the tree root */
                bitmap[tree[t]->yo][tree[t]->xo] = pBlack;
                bitmap[(unsigned) (tree[t]->yo - 1) % Sim->YMAX][tree[t]->xo] = pBlack;
                bitmap[(unsigned) (tree[t]->yo + 1) % Sim->YMAX][tree[t]->xo] = pBlack;
                bitmap[tree[t]->yo][(unsigned) (tree[t]->xo - 1) % Sim->XMAX] = pBlack;
                bitmap[tree[t]->yo][(unsigned) (tree[t]->xo + 1) % Sim->XMAX] = pBlack;
            }
            if (Sim->SHOWTREECENTER) { /* only a white dot at the center of growth */
                bitmap[tree[t]->y][tree[t]->x] = pWhite;
            }
        }
    }
}

void color_bitmap_ustory(void) {       
    /* This function applies to the UnderstoryCanopy_###.bmp files.
    */
    /* color the crown */
    int x, y, t;
    int num_trees = tree.size();
    for (x = 0; x < Sim->XMAX; x++)
        for (y = 0; y < Sim->YMAX; y++)
            /* canopy distribution */
            if (ustory_biomass[y][x] == -1) {
                bitmap[y][x] = pBlack;
            } else {
                bitmap[y][x] = tree[ustory_biomass[y][x]]->color;
            }
    /* black dot at the center of growth */
    for (t = 0; t < num_trees; t++)
        if (tree[t]->is_ustory) {
            bitmap[tree[t]->y][tree[t]->x].R = 0;
            bitmap[tree[t]->y][tree[t]->x].G = 0;
            bitmap[tree[t]->y][tree[t]->x].B = 0;
        }
}

void color_bitmap_ustory_hgt(void) {       
    int n, x, y, t;
    int num_trees = tree.size();
    short xx, yy;

    double gradient;

    /* color the crown */
    for (x = 0; x < Sim->XMAX; x++) {
        for (y = 0; y < Sim->YMAX; y++) {
            n = 0;
            if (ustory_biomass[y][x] == -1) {
                bitmap[y][x].R = 0x00;
                bitmap[y][x].G = 0x00;
                bitmap[y][x].B = 0x00;
            } else {
                /* height distribution */
                gradient = ((double)ustory_profile[y][x]) / (100.0 * ((double)Sim->HSCALE) * (1.0 + hRnd));
                bitmap[y][x].B =
                    (unsigned char) (tree[ustory_biomass[y][x]]->species->leaf.B * gradient);
                bitmap[y][x].G = 
                    (unsigned char) (tree[ustory_biomass[y][x]]->species->leaf.G * gradient);
                bitmap[y][x].R = 
                    (unsigned char) (tree[ustory_biomass[y][x]]->species->leaf.R * gradient);
                /* crown border points */
                for (yy = -1; yy <= 1; yy++) {
                    for (xx = -1; xx <= 1; xx++) {
                        if (ustory_biomass[(unsigned) (y + yy) % Sim->YMAX][(unsigned) (x + xx) % Sim->XMAX] == ustory_biomass[y][x]) {
                            n++;
                        }
                    }
                }
            }
            if (n > 3 && n < 8) {
                if (Sim->SHOWTREELIMITS) { /* crown border criterium */
                    bitmap[y][x] = pYellow; /* yellow line at the crown border */
                }
            }
        }
    }
    for (t = 0; t < num_trees; t++) {
        if (tree[t]->is_ustory) { /* tree is understory */
            if (Sim->SHOWTREEROOT) {  /* black cross at the tree root */
                bitmap[tree[t]->yo][tree[t]->xo] = pBlack;
                bitmap[(unsigned) (tree[t]->yo - 1) % Sim->YMAX][tree[t]->xo] = pBlack;
                bitmap[(unsigned) (tree[t]->yo + 1) % Sim->YMAX][tree[t]->xo] = pBlack;
                bitmap[tree[t]->yo][(unsigned) (tree[t]->xo - 1) % Sim->XMAX] = pBlack;
                bitmap[tree[t]->yo][(unsigned) (tree[t]->xo + 1) % Sim->XMAX] = pBlack;
            }
            if (Sim->SHOWTREECENTER) { /* white dot at the center of growth */
                bitmap[tree[t]->y][tree[t]->x] = pWhite;
            }
        }
    }
}

#ifdef __USE_ROOTS
void color_bitmap_rootstop(void) {       
    /* color the crown */
    int x, y, t, l, n, r = 0;
    int num_trees = tree.size();
    for (x = 0; x < Sim->XMAX; x++)
        for (y = 0; y < Sim->YMAX; y++) {
            for (l = 0, n = 0; l < ROOTMAX; l++)
                if (roots_top[y][x][l] != Sim->TMAX)
                    r++;
            /* root distribution */
            if (r == 1)
                bitmap[y][x] = tree[roots_top[y][x][0]].color;
            else if (r == 2)
                bitmap[y][x] = pLightGrey;
            else if (r == 3)
                bitmap[y][x] = pGrey;
            else if (r == 4)
                bitmap[y][x] = pDarkGrey;
            else if (r == 0)
                bitmap[y][x] = pBlack;
            r = 0;
            for (l = 0; l < ROOTMAX; l++)
                if (roots_top[y][x][l] != Sim->TMAX)
                    roots_top[y][x][l] = Sim->TMAX;
        }
    /* black dot at the center of growth */
    for (t = 0; t < num_trees; t++)
        if (tree[t]->d >= 0) {
            bitmap[roots[t].y][roots[t].x].R = 0xFF;
            bitmap[roots[t].y][roots[t].x].G = 0xFF;
            bitmap[roots[t].y][roots[t].x].B = 0xFF;
        }
}

void color_bitmap_rootsmid(void) {       
    /* color the crown */
    int x, y, t, l, r = 0;
    int num_trees = tree.size();
    for (x = 0; x < Sim->XMAX; x++)
        for (y = 0; y < Sim->YMAX; y++) {
            for (l = 0; l < ROOTMAX; l++)
                if (roots_mid[y][x][l] != Sim->TMAX)
                    r++;
            /* root distribution */
            if (r == 1)
                bitmap[y][x] = tree[roots_mid[y][x][0]].color;
            else if (r == 2)
                bitmap[y][x] = pLightGrey;
            else if (r == 3)
                bitmap[y][x] = pGrey;
            else if (r == 4)
                bitmap[y][x] = pDarkGrey;
            else if (r == 0)
                bitmap[y][x] = pBlack;
            r = 0;
            for (l = 0; l < ROOTMAX; l++)
                if (roots_mid[y][x][l] != Sim->TMAX)
                    roots_mid[y][x][l] = Sim->TMAX;
        }
    /* black dot at the center of growth */
    for (t = 0; t < num_trees; t++)
        if (tree[t]->d >= midqualifier) {
            bitmap[roots[t].y][roots[t].x].R = 0xFF;
            bitmap[roots[t].y][roots[t].x].G = 0xFF;
            bitmap[roots[t].y][roots[t].x].B = 0xFF;
        }
}

void color_bitmap_rootsbot(void) {       
    /* color the crown */
    int x, y, t, l, n, r = 0;
    int num_trees = tree.size();
    for (x = 0; x < Sim->XMAX; x++)
        for (y = 0; y < Sim->YMAX; y++) {
            for (l = 0, n = 0; l < ROOTMAX; l++)
                if (roots_bot[y][x][l] != Sim->TMAX)
                    r++;
            /* root distribution */
            if (r == 1)
                bitmap[y][x] = tree[roots_bot[y][x][0]].color;
            else if (r == 2)
                bitmap[y][x] = pLightGrey;
            else if (r == 3)
                bitmap[y][x] = pGrey;
            else if (r == 4)
                bitmap[y][x] = pDarkGrey;
            else if (r == 0)
                bitmap[y][x] = pBlack;
            r = 0;
            for (l = 0; l < ROOTMAX; l++)
                if (roots_bot[y][x][l] != Sim->TMAX)
                    roots_bot[y][x][l] = Sim->TMAX;
        }
    /* black dot at the center of growth */
    for (t = 0; t < num_trees; t++)
        if (tree[t]->d >= botqualifier) {
            bitmap[roots[t].y][roots[t].x].R = 0xFF;
            bitmap[roots[t].y][roots[t].x].G = 0xFF;
            bitmap[roots[t].y][roots[t].x].B = 0xFF;
        }
}
#endif

/* ================= */
/* SAVE PICTURE FILE */
/* ================= */
void save_png_crowns(int file_no) {
    int S,i;
    unsigned int t;
    double max_sec_val;
    int max_sec_i;
    double min_sec_val;
    int min_sec_i;
    /* create file name */
    std::stringstream ss;
    ss << std::setw(4) << file_no;
    std::string name = file_path + "Bmp/crowns_" + ss.str() + ".png";
    
    // compute max dimensions across species
    double max_height = 0;
    double max_width = 0;
    double max_trunk = 0;
    for (S = 0; S < Sim->SMAX; S++) {
        max_height = (max_height < Sim->species[S]->hMax) ? Sim->species[S]->hMax : max_height;
        max_width = (max_width < Sim->species[S]->wMax) ? Sim->species[S]->wMax : max_width;
        max_trunk = (max_trunk < Sim->species[S]->dMax) ? Sim->species[S]->dMax : max_trunk;
    }

    /* create new image */
    int indiv_width = 2.0*max_width + max_trunk+2;
    int total_width = indiv_width * 5;
    int indiv_height = max_height+1;
    int total_height = indiv_height * 3;
    pngwriter image(total_width, total_height, 0.0, name.c_str());

    //for (t = 0; t < tree.size(); t++) {
    for (unsigned int t1 = 0; t1 < 3; t1++) {
        for (unsigned int t2 = 0; t2 < 5; t2++) {
            t = t1*3+t2;
            int x_off = t2*indiv_width;
            int y_off = t1*indiv_height+1;
            if (t<tree.size() && tree[t]->is_ustory == 0) {
                // plot the max sector profile
                max_sec_i = std::distance(tree[t]->ws.begin(), std::max_element(tree[t]->ws.begin(), tree[t]->ws.end()));
                max_sec_val = tree[t]->ws[max_sec_i];
                for (i=0; i<max_sec_val-1; i++) {
                    image.line(
                            x_off+max_width+tree[t]->d/2.0+i, y_off+tree[t]->ComputeClaim((double) i, (double) max_sec_val,max_sec_i),
                            x_off+max_width+tree[t]->d/2.0+i+1, y_off+tree[t]->ComputeClaim((double) (i+1), (double) max_sec_val,max_sec_i),
                            ((int)(tree[t]->color.R)) / 255.0, 
                            ((int)(tree[t]->color.G)) / 255.0,
                            ((int)(tree[t]->color.B)) / 255.0);
                }
                // plot the min sector profile
                min_sec_i = std::distance(tree[t]->ws.begin(), std::min_element(tree[t]->ws.begin(), tree[t]->ws.end()));
                min_sec_val = tree[t]->ws[min_sec_i];
                for (i=0; i<min_sec_val-1; i++) {
                    image.line(
                            x_off+max_width-tree[t]->d/2.0-i, y_off+tree[t]->ComputeClaim((double) i, (double) min_sec_val,min_sec_i),
                            x_off+max_width-tree[t]->d/2.0-i-1, y_off+tree[t]->ComputeClaim((double) (i+1), (double) min_sec_val,min_sec_i),
                            ((int)(tree[t]->color.R)) / 255.0, 
                            ((int)(tree[t]->color.G)) / 255.0,
                            ((int)(tree[t]->color.B)) / 255.0);
                }
                // plot the trunk and finish the crown
                image.filledsquare(
                        x_off+max_width-tree[t]->d/2.0, y_off+0,
                        x_off+max_width+tree[t]->d/2.0, y_off+tree[t]->h,
                        ((int)(tree[t]->color.R)) / 255.0, 
                        ((int)(tree[t]->color.G)) / 255.0,
                        ((int)(tree[t]->color.B)) / 255.0);
            }
        }
    }

    char path[] = "include/cmtt10.ttf";
    char scale[] = "1m";
    image.plot_text(path, 12, 1, total_height - 15, 0, scale, 0.81,0.81,0.81);
//    image.line(5, total_height - 6, 15, total_height - 6,       0.81,0.81,0.81);
    image.line(5, total_height - 4, 15, total_height - 4,       0.81,0.81,0.81);

    image.close();
}


void save_png(std::string file_name, int file_no) {
    int x,y;
    /* create file name */
    std::stringstream ss;
    ss << std::setw(4) << file_no;
    std::string name = file_path + file_name + ss.str() + ".png";

    /* create new image */
    pngwriter image(Sim->XMAX, Sim->YMAX, 0.0, name.c_str());

    for (x = 0; x < Sim->XMAX; x++) {
        for (y = 0; y < Sim->YMAX; y++) {
            image.plot(x,y,((int)(bitmap[y][x].R)) / 255.0, 
                    ((int)(bitmap[y][x].G)) / 255.0,
                    ((int)(bitmap[y][x].B)) / 255.0);
        }
    }

    char path[] = "include/cmtt10.ttf";
    char scale[] = "1m";
    image.plot_text(path, 12, 10, 10, 0, scale, 1, 1, 1);
    image.line(15, 5, 25, 5,1,1,1);
    image.line(15, 6, 25, 6,1,1,1);
    image.line(15, 4, 25, 4,1,1,1);

    image.close();
}

void save_bmp(std::string file_name, int file_no) {
    FILE *stream;
    std::stringstream ss;

    /* create file name */
    ss << std::setw(4) << file_no;
    std::string name = file_path + file_name + ss.str() + ".bmp";

    if ((stream = fopen(name.c_str(), "wb+")) == NULL)
        std::cerr << "Error opening file " << name << " --> IGNORING\n";
    else {
        fputs("BM", stream);    /* Windows bitmap  */
        fwrite(&header, sizeof(struct DIB), 1, stream);
        //fwrite(bitmap, sizeof(struct PIXEL), Sim->XMAX * Sim->YMAX, stream);
        // CPP TODO FUCK
        fclose(stream);
    }
}




/* ============================= */
/* SAVE 3D MODEL (using POV-RAY) */
/* ============================= */
void save_3d(int file_no) {
    int sector;
    unsigned int t;
    /* create file connection */
    std::stringstream ss;
    ss << std::setw(4) << file_no;
    std::string name = file_path + "Bmp/scene_" + ss.str() + ".pov";
    std::ofstream stream;
    stream.open(name, std::ofstream::out);
    if (!stream.is_open()) {
        std::cerr << "Error opening file " << name << " --> IGNORING" << std::endl;
    } else {
        // camera + lights
        stream << "#include \"colors.inc\"\n\ncamera {\n				location <" << -0.5 * ((double) Sim->XMAX) << "," << -0.5 * ((double) Sim->YMAX)<< "," << ((double)Sim->XMAX + (double)Sim->YMAX) * 1.75 << ">\n				look_at <" << 0.5 * ((double) Sim->XMAX) << "," << 0.5 * ((double) Sim->YMAX)<< "," << ((double)Sim->XMAX + (double)Sim->YMAX) * 0.1875 << ">\n				up z\n				sky <0,0,1>\n}\n\nlight_source {\n				<" << 0.5 * ((double) Sim->XMAX) << "," << 0.5 * ((double) Sim->YMAX) << ",1000>\n				color rgb <1,1,1>\n				parallel\n				point_at <" << 0.5 * ((double) Sim->XMAX) << "," << 0.5 * ((double) Sim->YMAX) << ",0>\n}\n";
        // pigments
        stream << "#declare Pigment_1 = \npigment{ crackle turbulence 0.35 scale 0.45 \n         color_map{\n          [0.00 color rgb<0.25,0.5,0>]\n          [0.08 color rgb<0.75,1,0.5>]\n          [0.40 color rgb<0.5,1,0>]\n          [1.00 color rgb<0.75,1,0.5>]\n}\n}";
        stream << "#declare Pigment_2 = \npigment{ crackle turbulence 0.35 scale 0.45 \n         color_map{\n          [0.00 color rgb<1,1,0>]\n          [0.08 color rgb<1,1,0.5>]\n          [0.40 color rgb<1,1,0>] // true value\n          [1.00 color rgb<1,1,0.5>]\n}\n}";

            for (t=0; t<tree.size(); t++) {
                // output the trunks:
                double z = 1;
                stream << "cylinder {\n   <" << tree[t]->xo << "," << tree[t]->yo << "," << z << ">, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->h*0.8 << ">, " << tree[t]->d << "\n   pigment { color rgb <0.3,0.2,0> }\n}\n";
                // output the crowns:
                std::string color;
                if (tree[t]->species->shadetolerance == 0) {
                    color = "texture{ pigment{ Pigment_1  }           finish { phong 1 }         }";
                } else {
                    color = "texture{ pigment{ Pigment_2  }           finish { phong 1 }         }";
                }
                //				color="<1,0,0>";
                for (sector=0; sector<tree[t]->SECTORS; sector++) {
                    double x1 = tree[t]->ws[sector] * cos(M_PI * (sector/4.0)) + tree[t]->x;
                    double y1 = tree[t]->ws[sector] * sin(M_PI * (sector/4.0)) + tree[t]->y;
                    double x2 = tree[t]->ws[sector] * cos(M_PI * ((sector+1)/4.0)) + tree[t]->x;
                    double y2 = tree[t]->ws[sector] * sin(M_PI * ((sector+1)/4.0)) + tree[t]->y;
                    // crown
                    stream << "polygon {\n4, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->h << ">, <" << x1 << "," << y1 << "," << z+tree[t]->cbhs[sector] << ">, <" << x2 << "," << y2 << "," << z+tree[t]->cbhs[sector] << ">, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->h << "> \n" << color << " \n}\n";
                    // vert wall 1
                    stream << "polygon {\n4, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->h << ">, <" << x1 << "," << y1 << "," << z+tree[t]->cbhs[sector] << ">, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->cbhs[sector] << ">, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->h << "> \n" << color << " \n}\n";
                    // vert wall 2
                    stream << "polygon {\n4, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->h << ">, <" << x2 << "," << y2 << "," << z+tree[t]->cbhs[sector] << ">, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->cbhs[sector] << ">, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->h << "> \n" << color << " \n}\n";
                    // floor:
                    stream << "polygon {\n4, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->cbhs[sector] << ">, <" << x1 << "," << y1 << "," << z+tree[t]->cbhs[sector] << ">, <" << x2 << "," << y2 << "," << z+tree[t]->cbhs[sector] << ">, <" << tree[t]->x << "," << tree[t]->y << "," << z+tree[t]->cbhs[sector] << "> \n" << color << " \n}\n";
                    //              stream << "bicubic_patch {\n    type 0\n    flatness 0.01\n   u_steps 4\n   v_steps 4\n     ";
                    //              //first wall is the bottom:
                    //              stream << "<0, 0, 2>, <1, 0, 0>, <2, 0, 0>, <3, 0,-2>,\n     <0, 1  0>, <1, 1, 0>, <2, 1, 0>, <3, 1, 0>,\n      <0, 2, 0>, <1, 2, 0>, <2, 2, 0>, <3, 2, 0>,\n     <0, 3, 2>, <1, 3, 0>, <2, 3, 0>, <3, 3, -2>\n}";
                }
            }
        stream.close();
    }
}



/* ============= */
/* SAVE RAW DATA */
/* ============= */
void save_data(std::string file_name, int file_no) {
    //FILE *stream;
    //unsigned int num_trees = tree.size();
    //unsigned int i;

    // CPP
    std::cout << "FIXME: no save_data\n";
    //
    //		/* create file name */
    //        std::stringstream ss;
    //        ss << std::setw(4) << file_no;
    //		std::string name = file_path + file_name + ss.str() + ".dat";
    //
    //		if ((stream = fopen(name.c_str(), "wb")) == NULL)
    //			std::cerr << "Error opening file " << name << " --> IGNORING\n";
    //		else {
    //			fwrite(&num_trees, sizeof(num_trees), 1, stream);
    //			/* fwrite(tree, sizeof(struct UNIT), num_trees, stream); */
    //			for (i = 0; i < num_trees; i++)
    //				fwrite((tree + i), sizeof(struct UNIT), 1, stream);
    //			fclose(stream);
    //		}
}

/* ====================== */
/* BASIC CANOPY TREE DATA */
/* ====================== */
void save_alltrees(int i) {
    int s;
    unsigned int t;
    for (t = 0; t < tree.size(); t++) {
        file_all_trees << i << "\t" // year
            << tree[t]->n << "\t"   /* tree index */
            << tree[t]->x << "\t"  /* jean: x and y positions of the crown, in dm */
            << tree[t]->y << "\t"
            << tree[t]->s << "\t"  /* tree species index */ // TODO: handle the intermediate species case
            << tree[t]->d * 10.0 << "\t" /* trunk diameter in cm */
            << tree[t]->h / 10.0 << "\t" /* tree height in m */
            << (double) tree[t]->ca / 100.0 << "\t"  /* canopy area in m2 */
            << (double) tree[t]->water_uptake / tree[t]->water_needed << "\t"  /* water ratio */
            << tree[t]->xo << "\t"  /* jean: x and y positions of the tree roots (i.e. of the seed), in dm */
            << tree[t]->yo << "\t";
        for (s = 0; s < Sim->SECTORS; s++) {
            file_all_trees << (double) tree[t]->cbhs[s] / 10.0 << "\t"; /* sectors' cbh */
        }
        for (s = 0; s < Sim->SECTORS; s++) {
            file_all_trees << (double) tree[t]->ws[s] / 10.0 << "\t"; /* sectors' cbh */
        }
        file_all_trees << std::endl;
    }
    file_all_trees.flush();
}

void save_stats(int i) {
    unsigned int t;
    int x, y;
    double basalarea=0.0;
    double shadetoleranceindex_nb = 0.0;
    double shadetoleranceindex_ba = 0.0;
    double totalbasalarea = 0.0;
    double age = 0.0;
    double meanassymetry = 0.0, tmp_a, tmp_A, tmp_ws1, tmp_ws2, tmp;
    double occupied_canopy = 0.0;
    int realtreenb = 0;
    int understorynb = 0;
    int overstorynb = 0;
    double understoryba = 0;
    double overstoryba = 0;
    double understoryage = 0;
    double overstoryage = 0;
    // Added in Nov. 2016: water uptake/available ratio ("WATER_CONSUMED")
    // and water needed/available ratio ("WATER_NEEDED")
    double available_water = Sim->XMAX * Sim->YMAX * Sim->WATER_AVAILABLE;
    double sucked_water = 0.;
    double needed_water = 0.;
    double sucked_water_a = 0.;
    double sucked_water_b = 0.;
    double sucked_water_c = 0.;
    double sucked_water_d = 0.;
    double sucked_water_e = 0.;
    double needed_water_a = 0.;
    double needed_water_b = 0.;
    double needed_water_c = 0.;
    double needed_water_d = 0.;
    double needed_water_e = 0.;
    int pop_a = 0;
    int pop_b = 0;
    int pop_c = 0;
    int pop_d = 0;
    int pop_e = 0;
    int pop_intolerant = 0;
    double sucked_water_tolerant = 0.;
    double sucked_water_intolerant = 0.;
    double needed_water_tolerant = 0.;
    double needed_water_intolerant = 0.;
    int pop_tolerant = 0;

    for (t = 0; t < tree.size(); t++) {
        tmp_a = -1;
        tmp_A = -1;
        sucked_water += tree[t]->water_uptake;
        needed_water += tree[t]->water_needed;
        if (tree[t]->d < 0.5) {
            sucked_water_a += tree[t]->water_uptake;
            needed_water_a += tree[t]->water_needed;
            pop_a++;
        } else if (tree[t]->d < 1.) {
            sucked_water_b += tree[t]->water_uptake;
            needed_water_b += tree[t]->water_needed;
            pop_b++;
        } else if (tree[t]->d < 1.5) {
            sucked_water_c += tree[t]->water_uptake;
            needed_water_c += tree[t]->water_needed;
            pop_c++;
        } else if (tree[t]->d < 2.) {
            sucked_water_d += tree[t]->water_uptake;
            needed_water_d += tree[t]->water_needed;
            pop_d++;
        } else {
            sucked_water_e += tree[t]->water_uptake;
            needed_water_e += tree[t]->water_needed;
            pop_e++;
        }
        if (tree[t]->species->shadetolerance == 1) {
            // shade tolerant
            sucked_water_tolerant += tree[t]->water_uptake;
            needed_water_tolerant += tree[t]->water_needed;
            pop_tolerant++;
        } else {
            // shade intolerant
            sucked_water_intolerant += tree[t]->water_uptake;
            needed_water_intolerant += tree[t]->water_needed;
            pop_intolerant++;
        }
        for (int sector=0; sector < tree[t]->SECTORS/2; sector++) {
            tmp_ws1 = ((tree[t]->ws[sector] > tree[t]->ustory_ws[sector]) ? tree[t]->ws[sector]: tree[t]->ustory_ws[sector]);
            tmp_ws2 = ((tree[t]->ws[sector+tree[t]->SECTORS/2-1] > tree[t]->ustory_ws[sector+tree[t]->SECTORS/2-1]) ? tree[t]->ws[sector+tree[t]->SECTORS/2-1]: tree[t]->ustory_ws[sector+tree[t]->SECTORS/2-1]);
            tmp = tmp_ws1 + tmp_ws2;
            if (tmp_a > tmp || tmp_a == -1) {
                tmp_a = tmp;
            }
            if (tmp_A < tmp || tmp_A == -1) {
                tmp_A = tmp;
            }
        }
        if (tree[t]->d>0.7 && tree[t]->is_ustory != 1) {
            // we exclude trees that are too small or that have no crown at all
            basalarea = M_PI * (tree[t]->d * tree[t]->d) / 4.0;
            totalbasalarea += basalarea;
            shadetoleranceindex_nb += tree[t]->species->shadetolerance;
            shadetoleranceindex_ba += tree[t]->species->shadetolerance * basalarea;
            age += tree[t]->age;
            meanassymetry += tmp_a/tmp_A;
            realtreenb++;
        }
        if (tree[t]->is_ustory == 1) {
            understorynb++;
            understoryba += M_PI * (tree[t]->d * tree[t]->d) / 4.0;
            understoryage += tree[t]->age;
        } else {
            overstorynb++;
            overstoryba += M_PI * (tree[t]->d * tree[t]->d) / 4.0;
            overstoryage += tree[t]->age;
        }
    }

    for (x=0; x < Sim->XMAX; x++) {
        for (y=0; y < Sim->YMAX; y++) {
            if (biomass[y][y] >= 0) {
                occupied_canopy++;
            }
        }
    }
    double meanage;
    if (age != 0.) {
        meanage = age / ((double)realtreenb);
    } else {
        meanage = 0.;
    }
    occupied_canopy = occupied_canopy / ((double)Sim->XMAX) / ((double)Sim->YMAX);

    file_stats << i << "\t" // year
            << realtreenb << "\t"   /* total overstory pop */
            << meanage << "\t"   /* mean age */
            << totalbasalarea << "\t"  /* basal area per unit */
            << shadetoleranceindex_nb / (realtreenb) << "\t" /* STI based on tree number */
            << shadetoleranceindex_ba / totalbasalarea << "\t" /* STI based on tree number */
            << meanassymetry / ((double)realtreenb) << "\t" /* mean assymetry */
            << occupied_canopy << "\t" /* occupied canopy ratio */
            << understorynb << "\t" /* total understory pop */
            << overstorynb << "\t" /* total overstory pop */
            << understoryba << "\t" /* total understory ba */
            << overstoryba << "\t" /* total overstory ba */
            << understoryage / ((double)understorynb) << "\t" /* total understory age */
            << overstoryage / ((double)overstorynb) << "\t" /* total overstory age */
            << sucked_water / available_water << "\t" /* consumed water ratio */
            << needed_water / available_water << "\t" /* water fulfillment ratio */
            << sucked_water_a / available_water<< "\t"
            << sucked_water_b / available_water<< "\t"
            << sucked_water_c / available_water<< "\t"
            << sucked_water_d / available_water<< "\t"
            << sucked_water_e / available_water<< "\t"
            << needed_water_a / available_water<< "\t"
            << needed_water_b / available_water<< "\t"
            << needed_water_c / available_water<< "\t"
            << needed_water_d / available_water<< "\t"
            << needed_water_e / available_water<< "\t"
            << sucked_water_intolerant / available_water<< "\t"
            << needed_water_intolerant / available_water<< "\t"
            << sucked_water_tolerant / available_water<< "\t"
            << needed_water_tolerant / available_water<< "\t"
            << pop_a << "\t"
            << pop_b << "\t"
            << pop_c << "\t"
            << pop_d << "\t"
            << pop_d << "\t"
            << pop_e << "\t"
            << pop_tolerant << "\t"
            << pop_intolerant
            << std::endl;
    file_stats.flush();
    std::cout << i << "\t"
            << tree.size() << "\t"
            << meanage << "\t"
            << totalbasalarea << "\t"
            << shadetoleranceindex_nb / ((double)realtreenb) << "\t"
            << shadetoleranceindex_ba / totalbasalarea << "\t"
            << meanassymetry / ((double)realtreenb) << "\t"
            << occupied_canopy << "\t"
            << understorynb << "\t"
            << overstorynb << "\t"
            << understoryba << "\t"
            << overstoryba << "\t"
            << understoryage / ((double)understorynb) << "\t"
            << overstoryage / ((double)overstorynb) << "\t"
            << sucked_water / available_water << "\t"
            << needed_water / available_water
            << std::endl;
}


void save_tree(std::string file_name, int file_no) {
    FILE *stream;
    int t;
    int s;
    int num_trees = tree.size();

    std::ofstream outfile ("test.txt");

    /* create file name */
    std::stringstream ss;
    ss << std::setw(4) << file_no;
    std::string name = file_path + file_name + ss.str() + ".txt";

    if ((stream = fopen(name.c_str(), "w+")) == NULL)
        std::cerr << "Error opening file " << name << " --> IGNORING\n";
    else {
        //fprintf(stream, "ID \n");
        for (t = 0; t < num_trees; t++) {
            //				if (tree[t]->ca > 0) {
            //fprintf(stream, "%4i\t%4i\t%4i\t%4i\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%4i\t%4i\n", tree[t]->n,   /* tree index */
            //		tree[t]->x,  /* jean: x and y positions of the crown, in dm */
            //		tree[t]->y,
            //		tree[t]->s,  /* tree species index */
            //		tree[t]->d * 10, /* trunk diameter in cm */
            //		tree[t]->h / 10, /* tree height in m */
            //		(double) tree[t]->ca / 100,  /* canopy area in m2 */
            //		(double) tree[t]->water_uptake / tree[t]->water_needed,  /* water ratio */
            //		(double) tree[t]->cbh / 10,  /* crown base height in m */
            //		(double) tree[t]->brdl / 10, /* crown border length in m */
            //		tree[t]->brdh / 10,  /* crown border average height in m */
            //		cross_section(*tree[t], hOver) / 100,
            //		tree[t]->xo,  /* jean: x and y positions of the tree roots (i.e. of the seed), in dm */
            //		tree[t]->yo
            //		);   /* crown cross-section at overlap height in m2 */
            for (s = 0; s < Sim->SECTORS; s++)
                fprintf(stream, "%6.2f\t", (double) tree[t]->cbhs[s] / 10); /* sectors' cbh */
            for (s = 0; s < Sim->SECTORS; s++)
                fprintf(stream, "%6.2f\t", (double) tree[t]->ws[s]);
            fprintf(stream, "\n");
            //				}
        }
        fclose(stream);
    }
}

/* ================ */
/* BASIC MODEL DATA */
/* ================ */
//void save_balance(std::string file_name, int file_no) {
//    FILE *stream;
//    std::string name = file_path + file_name + ".txt";
//    int s;
//    if (file_no == 1) {
//        stream = fopen(name.c_str(), "w");
//    } else {
//        stream = fopen(name.c_str(), "a+");
//    }
//    if (stream == NULL) {
//        std::cerr << "Error opening file " << name << " --> IGNORING\n";
//    } else {
//        for (s = 0; s <= Sim->SMAX; s++)
//            fprintf(stream, "%4i\t%4i\t+%4i\t-%4i\t=%4i\t(%4i)\t%5.2f\n", file_no, s,   /* spesies No */
//                    stats[s].tNew,   /* number of trees planted last cycle */
//                    stats[s].tDead,  /* number of trees died last cycle */
//                    stats[s].tNow,   /* total number of trees */
//                    stats[s].tActive,    /* number of trees with a canopy */
//                    hOver / 10);    /* projected crown cross-section overlap height in m (approximator) */
//        fclose(stream);
//    }
//}

/* ================ */
/* 3D SCENE EXPORTS */
/* ================ */
void save_space(std::string file_name, int file_no) {
    FILE *stream;
    short x,y;
    //double plan_alt_max = ((double)Sim->XMAX) * sin(ANGLE);

    /* create file name */
    std::stringstream ss;
    ss << std::setw(4) << file_no;
    std::string name = file_path + file_name + ss.str() + ".txt";

    if ((stream = fopen(name.c_str(), "w+")) == NULL)
        std::cerr << "Error opening file " << name << " --> IGNORING\n";
    else {
        for (x = 0; x < Sim->XMAX; x++) {
            for (y = 0; y < Sim->YMAX; y++) {
                //					if (xyclaims[x][y]->z <= 0) {
                //							xyclaims[x][y]->z = 0;
                //							xyclaims[x][y]->is_sunny = 0;
                //					}// else {
                fprintf(stream, "%4i\t%4i\t%f\t%.2f\t%.2f\t%.2f\n",
                        x,
                        y,
                        xyclaims[x][y]->z,
                        //								tree[xyclaims[x][y]->t].s // to output the specie
                        //1.0 / (1.0 + (double) (xyclaims[x][y]->is_sunny))
                        (double)xyclaims[x][y]->is_sunny,
                        get_altitude((double)x,(double)y),
                        xyclaims[x][y]->v
                       );   
                //}
            }
            fprintf(stream, "\n");
        }
        fclose(stream);
    }
}

#if VORONOI
/* If we want to generate the output voronoi diagrams, set VORONOI to 1 in mod.h
 * This will generate text files that will be used in voronoi.cpp found in
 * Work/Models/voronoi */
void voronoi(int y) {
    FILE *wstream = NULL;
    FILE *cstream = NULL;
    std::string buffer[10];
    std::string wfile_name[1000];
    std::string cfile_name[1000];
    std::string *w = (char *)"weights_";
    std::string *c = (char *)"positions_";
    std::string *flag =(char *) "w+";
    int i = 1;
    const int ntrees = stats[Sim->SMAX].tNow + stats[Sim->SMAX].ustNow;

    std::strcpy(wfile_name, "Work/Models/voronoi/");
    std::strcpy(cfile_name, "Work/Models/voronoi/");
    std::strcat(wfile_name, w);
    std::strcat(cfile_name, c);
    std::sprintf(buffer, "%04i", y);
    std::strcat(wfile_name, buffer);
    std::strcat(cfile_name, buffer);
    std::strcat(wfile_name, ".txt");
    std::strcat(cfile_name, ".txt");
    if (y == 1)
        flag =(char *) "w";
    if ((NULL != (wstream = fopen(wfile_name, flag))) &&
            (NULL != (cstream = fopen(cfile_name, flag)))) {
        fprintf(cstream, "%d\n", ntrees);
        for (; i < ntrees; i++) {
            if (tree[i]->is_ustory == 0) {
                fprintf(cstream, "%d,%d,\n", tree[i]->xo, tree[i]->yo);
                fprintf(wstream, "%6.2f\n", tree[i]->d);
            }
        }
        fclose(wstream);
        fclose(cstream);
    }
}
#endif

void cancel_key(int signo) {
    freeall();
    printf("\n");
    exit(0);
}

/* =================== */
/* OBSOLETE PROCEDURES */
/* =================== */
/*void save_heightarea(char *txt_name, int file_no) {
 *    FILE *stream;
 *
 *    char buffer[3];
 *
 *    char file_name[100];
 *
 *    strcpy(file_name, file_path);
 *    strcat(file_name, txt_name);
 *    sprintf(buffer, "%03i", file_no);
 *    strcat(file_name, buffer);
 *    strcat(file_name, ".txt");
 *    if ((stream = fopen(file_name, "w")) == NULL)
 *        printf("Error opening file %s\n", file_name);
 *    else {
 *        for (unsigned short t = 0; t < stats[Sim->SMAX].tNow; t++)
 *            if (tree[t]->ca > 0)
 *                fprintf(stream, "%4i\t%6.2f\t%8.2f\t%i\n", tree[t]->n,
 *                        tree[t]->d, tree[t]->h, tree[t]->ca);
 *        fclose(stream);
 *    }
 *}
 *void save_profile(char *file_name, int file_no) {
 *    FILE *stream;
 *
 *    char buffer[3];
 *
 *    char name[100];
 *
 *    strcpy(name, file_path);
 *    strcat(name, file_name);
 *    sprintf(buffer, "%03i", file_no);
 *    strcat(name, buffer);
 *    strcat(name, ".txt");
 *    if ((stream = fopen(name, "w")) == NULL)
 *        printf("Error opening file %s\n", name);
 *    else {
 *        for (unsigned short x = 0; x < xMax; x++)
 *            for (unsigned short y = 0; y < yMax; y++)
 *                fprintf(stream, "%8.3f\n", (double) profile[y][x] / 1000);
 *    }
 *    fclose(stream);
 *}
 */



