#include "simulation.hpp"
#include "species.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iterator>
#include <fstream>
#include <sstream>
#include <iostream>

#ifndef GIT_VERSION
#include "gitversion_tmpfile.h"
#endif

double Simulation::initialize_nSeedTotal()
{
  double ratio = ((double)XMAX*(double)YMAX)/(1024.0*1024.0);
  return (SEEDS_PER_HA * ratio);
}


Simulation::Simulation(int ac, char** av) {
  try {
    std::string opt;
    std::string config_file;

    // Declare a group of options that will be 
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
      ("version,v", "print version string")
      ("help", "produce help message")
      ("config,c", po::value<std::string>(&config_file)->default_value("default.ini"),
       "name of a file of a configuration.")
      ;

    // Declare a group of options that will be 
    // allowed both on command line and in
    // config file
    po::options_description config("Output directory");
    config.add_options()
      ("output,o", po::value<std::string>(&opt)->default_value("default"), 
       "base directory for output")
      ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
      // SIMULATION OPTIONS
      ("HSCALE", po::value< int >()->default_value(250), "HSCALE")
      ("CYCLES", po::value< int >()->default_value(1500), "CYCLES")
      ("XMAX", po::value< int >()->default_value(512), "XMAX")
      ("YMAX", po::value< int >()->default_value(512), "YMAX")
      ("ZMAX", po::value< int >()->default_value(512), "ZMAX")
      ("WATER_AVAILABLE", po::value< double >()->default_value(1.0), "WATER_AVAILABLE")
      ("SEED", po::value< int >()->default_value(42), "SEED")
      ("SECTORS", po::value< int >()->default_value(8), "SECTORS")
      ("CORIENTATION", po::value< double >()->default_value(0.5*M_PI), "CORIENTATION")
      ("CANGLE", po::value< double >()->default_value(0), "CANGLE")
      ("CSUNALT", po::value< double >()->default_value(0.5*M_PI), "CSUNALT")
      ("CSUNORIENTATION", po::value< double >()->default_value(1.54*M_PI), "CSUNORIENTATION")
      ("SMAX", po::value< int >()->default_value(2), "SMAX")
      ("SEEDS_PER_HA", po::value< double >()->default_value(40.0), "SEEDS_PER_HA")
      ("COMPLETEDISTURBANCE", po::value< double >()->default_value(0.0), "COMPLETEDISTURBANCE") // probability of complete disturbance (default: 0, no complete disturbance)
      ("ROOT", po::value< int >()->default_value(0), "ROOT") // default: 0 (no root)
      ("LIGHTWATERGROWTH", po::value< int >()->default_value(1), "LIGHTWATERGROWTH") // default: 1 (multiplicative influence)
      ("LIGHTWATERMORTALITY", po::value< int >()->default_value(3), "LIGHTWATERMORTALITY") // default: 3 (additive influence)
      ("RESTRICTCROWN", po::value< int >()->default_value(0), "RESTRICTCROWN")
      ("VERTICALGROWTH", po::value< int >()->default_value(0), "VERTICALGROWTH") // when set to 0, this neglects hinc
      ("RESTRICTHEIGHT", po::value< int >()->default_value(2), "RESTRICTHEIGHT")
      ("RESTRICTRADIUS", po::value< int >()->default_value(0), "RESTRICTRADIUS")
      ("POTENTIALCROWN", po::value< int >()->default_value(1), "POTENTIALCROWN")
//      ("ZRnd", po::value< double >()->default_value(1), "ZRnd")
      ("INTERMEDIATESPECIES", po::value< int >()->default_value(0), "INTERMEDIATESPECIES")
      ("SEEDANYWHERE", po::value< int >()->default_value(1), "SEEDANYWHERE")
      ("CROWNSHAPE", po::value< int >()->default_value(1), "CROWNSHAPE")
      ("REAL3D", po::value< int >()->default_value(0), "REAL3D")
      ("OUTPUTMOD", po::value< int >()->default_value(1), "OUTPUTMOD")
      ("OUTPUTBMP", po::value< int >()->default_value(1), "OUTPUTBMP")
      ("OUTPUTALLTREES", po::value< int >()->default_value(1), "OUTPUTALLTREES")
      ("OUTPUTSTATS", po::value< int >()->default_value(1), "OUTPUTSTATS")
      ("SHOWTREELIMITS", po::value< int >()->default_value(0), "SHOWTREELIMITS")
      ("SHOWTREEROOT", po::value< int >()->default_value(1), "SHOWTREEROOT")
      ("SHOWTREECENTER", po::value< int >()->default_value(1), "SHOWTREECENTER")
      ("REDO_XP", po::value< int >()->default_value(0), "REDO_XP") // default: 0 (ask)
      // S0 is typically a shade tolerant species
      ("S0.mrShade", po::value< double >(), "S0.mrShade")
      ("S0.rinc", po::value< double >(), "S0.rinc")
      ("S0.hinc", po::value< double >(), "S0.hinc")
      ("S0.growth_factor", po::value< double >(), "S0.growth_factor")
      // S1 is typically a shade intolerant species
      ("S1.mrShade", po::value< double >(), "S1.mrShade")
      ("S1.rinc", po::value< double >(), "S1.rinc")
      ("S1.hinc", po::value< double >(), "S1.hinc")
      ("S1.growth_factor", po::value< double >(), "S1.growth_factor")
      // most of these have already default values specified in Species' constructor, but we still must to declare here those that we want to modify
      ("S0.mrLight", po::value< double >(), "S0.mrLight")
      ("S1.mrLight", po::value< double >(), "S1.mrLight")
      ("S2.mrLight", po::value< double >(), "S2.mrLight")
      ("S0.water_mr", po::value< double >(), "S0.water_mr")
      ("S1.water_mr", po::value< double >(), "S1.water_mr")
      ("S2.water_mr", po::value< double >(), "S2.water_mr")
      ("S0.water_needed", po::value< double >(), "S0.water_needed")
      ("S1.water_needed", po::value< double >(), "S1.water_needed")
      ("S2.water_needed", po::value< double >(), "S2.water_needed")
      ("S0.leafR", po::value< int >(), "S0.leafR")
      ("S0.leafG", po::value< int >(), "S0.leafG")
      ("S0.leafB", po::value< int >(), "S0.leafB")
      ("S1.leafR", po::value< int >(), "S1.leafR")
      ("S1.leafG", po::value< int >(), "S1.leafG")
      ("S1.leafB", po::value< int >(), "S1.leafB")
      ("S2.leafR", po::value< int >(), "S2.leafR")
      ("S2.leafG", po::value< int >(), "S2.leafG")
      ("S2.leafB", po::value< int >(), "S2.leafB")
      ("S0.shadetolerance", po::value< double >(), "S0.shadetolerance")
      ("S1.shadetolerance", po::value< double >(), "S1.shadetolerance")
      ("S2.shadetolerance", po::value< double >(), "S2.shadetolerance")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);

    po::positional_options_description p;
    p.add("input-file", -1);

    po::variables_map vm;
    store(po::command_line_parser(ac, av).
        options(cmdline_options).positional(p).run(), vm);
    notify(vm);

    std::ifstream ifs(config_file.c_str());
    if (!ifs)
    {
      std::cout << std::endl << "!!!!!!!!!" << std::endl << "WARNING: can not open config file: " << config_file << std::endl << "!!!!!!!!!!!!!" << std::endl;
    }
    else
    {
      store(parse_config_file(ifs, config_file_options), vm);
      notify(vm);
    }

    if (vm.count("help")) {
      std::cout << visible << "\n";
    }

    if (vm.count("version")) {
      std::cout << "Git commit: " << GIT_VERSION << std::endl;
    }

    if (vm.count("HSCALE"))
    {
      std::cout << "HSCALE is: " << vm["HSCALE"].as< int >() << "\n";
      HSCALE = vm["HSCALE"].as< int >();
    }
    //if (vm.count("HSCALE")) {HSCALE = vm["HSCALE"].as< int >();}
    if (vm.count("CYCLES")) {CYCLES = vm["CYCLES"].as< int >();}
    if (vm.count("XMAX")) {XMAX = vm["XMAX"].as< int >();}
    if (vm.count("YMAX")) {YMAX = vm["YMAX"].as< int >();}
    if (vm.count("ZMAX")) {ZMAX = vm["ZMAX"].as< int >();}
    if (vm.count("WATER_AVAILABLE")) {WATER_AVAILABLE = vm["WATER_AVAILABLE"].as< double >();}
    if (vm.count("SEED")) {SEED = vm["SEED"].as< int >();}
    if (vm.count("SECTORS")) {SECTORS = vm["SECTORS"].as< int >();}
    if (vm.count("CORIENTATION")) {CORIENTATION = vm["CORIENTATION"].as< double >();}
    if (vm.count("CANGLE")) {CANGLE = vm["CANGLE"].as< double >();}
    if (vm.count("CSUNALT")) {CSUNALT = vm["CSUNALT"].as< double >();}
    if (vm.count("CSUNORIENTATION")) {CSUNORIENTATION = vm["CSUNORIENTATION"].as< double >();}
    if (vm.count("SMAX")) {SMAX = vm["SMAX"].as< int >();}
    if (vm.count("SEEDS_PER_HA")) {SEEDS_PER_HA = vm["SEEDS_PER_HA"].as< double >();}
    if (vm.count("COMPLETEDISTURBANCE")) {COMPLETEDISTURBANCE = vm["COMPLETEDISTURBANCE"].as< double >();}
    if (vm.count("ROOT")) {ROOT = vm["ROOT"].as< int >();}
    if (vm.count("LIGHTWATERGROWTH")) {LIGHTWATERGROWTH = vm["LIGHTWATERGROWTH"].as< int >();}
    if (vm.count("LIGHTWATERMORTALITY")) {LIGHTWATERMORTALITY = vm["LIGHTWATERMORTALITY"].as< int >();}
    if (vm.count("RESTRICTCROWN")) {RESTRICTCROWN = vm["RESTRICTCROWN"].as< int >();}
    if (vm.count("VERTICALGROWTH")) {VERTICALGROWTH = vm["VERTICALGROWTH"].as< int >();}
    if (vm.count("RESTRICTHEIGHT")) {RESTRICTHEIGHT = vm["RESTRICTHEIGHT"].as< int >();}
    if (vm.count("RESTRICTRADIUS")) {RESTRICTRADIUS = vm["RESTRICTRADIUS"].as< int >();}
    if (vm.count("POTENTIALCROWN")) {POTENTIALCROWN = vm["POTENTIALCROWN"].as< int >();}
//    if (vm.count("ZRnd")) {ZRnd = vm["ZRnd"].as< double >();}
    if (vm.count("INTERMEDIATESPECIES")) {INTERMEDIATESPECIES = vm["INTERMEDIATESPECIES"].as< int >();}
    if (vm.count("SEEDANYWHERE")) {SEEDANYWHERE = vm["SEEDANYWHERE"].as< int >();}
    if (vm.count("CROWNSHAPE")) {CROWNSHAPE = vm["CROWNSHAPE"].as< int >();}
    if (vm.count("REAL3D")) {REAL3D = vm["REAL3D"].as< int >();}
    if (vm.count("OUTPUTMOD")) {OUTPUTMOD = vm["OUTPUTMOD"].as< int >();}
    if (vm.count("OUTPUTALLTREES")) {OUTPUTALLTREES = vm["OUTPUTALLTREES"].as< int >();}
    if (vm.count("OUTPUTSTATS")) {OUTPUTSTATS = vm["OUTPUTSTATS"].as< int >();}
    if (vm.count("OUTPUTBMP")) {OUTPUTBMP = vm["OUTPUTBMP"].as< int >();}
    if (vm.count("SHOWTREELIMITS")) {SHOWTREELIMITS = vm["SHOWTREELIMITS"].as< int >();}
    if (vm.count("SHOWTREEROOT")) {SHOWTREEROOT = vm["SHOWTREEROOT"].as< int >();}
    if (vm.count("SHOWTREECENTER")) {SHOWTREECENTER = vm["SHOWTREECENTER"].as< int >();}
    if (vm.count("REDO_XP")) {REDO_XP = vm["REDO_XP"].as< int >();}

    output_dirname.assign(opt);

    // SPECIES INITIALIZATION:
    for (int i=0; i<SMAX; i++) {
      std::cout << "Initializing species number " << i << std::endl;
      //species.push_back(Species(i)); 
      species.push_back(std::shared_ptr<Species> (new Species(i)));
      {
        std::stringstream ss;
        ss << "S" << i << ".mrShade";
        if (vm.count(ss.str().c_str()))
          species[i]->mrShade = vm[ss.str().c_str()].as< double >();
      }
      {
        std::stringstream ss;
        ss << "S" << i << ".mrLight";
        if (vm.count(ss.str().c_str()))
          species[i]->mrLight = vm[ss.str().c_str()].as< double >();
      }
      {
        std::stringstream ss;
        ss << "S" << i << ".rinc";
        if (vm.count(ss.str().c_str()))
          species[i]->rinc = vm[ss.str().c_str()].as< double >();
      }
      {
        std::stringstream ss;
        ss << "S" << i << ".hinc";
        if (vm.count(ss.str().c_str()))
          species[i]->hinc = vm[ss.str().c_str()].as< double >();
      }
      {
        std::stringstream ss;
        ss << "S" << i << ".growth_factor";
        if (vm.count(ss.str().c_str()))
          species[i]->growth_factor = vm[ss.str().c_str()].as< double >();
      }
      {
        std::stringstream ss;
        ss << "S" << i << ".leafR";
        if (vm.count(ss.str().c_str()))
          species[i]->leaf.R = (unsigned char) vm[ss.str().c_str()].as< int >();
      }
      {
        std::stringstream ss;
        ss << "S" << i << ".leafG";
        if (vm.count(ss.str().c_str()))
          species[i]->leaf.G = (unsigned char) vm[ss.str().c_str()].as< int >();
      }
      {
        std::stringstream ss;
        ss << "S" << i << ".leafB";
        if (vm.count(ss.str().c_str())) {
          species[i]->leaf.B = (unsigned char) vm[ss.str().c_str()].as< int >();
        }
      }
      {
        std::stringstream ss;
        ss << "S" << i << ".water_mr";
        if (vm.count(ss.str().c_str())) {
          species[i]->water_mr = vm[ss.str().c_str()].as< double >();
          std::cout << "Using the water mortality: " << species[i]->water_mr << std::endl;
        }
      }
      {
        std::stringstream ss;
        ss << "S" << i << ".water_needed";
        if (vm.count(ss.str().c_str())) {
          species[i]->water_needed = vm[ss.str().c_str()].as< double >();
          std::cout << "Using the water need: " << species[i]->water_needed << std::endl;
        }
      }
      {
        std::stringstream ss;
        ss << "S" << i << ".shadetolerance";
        if (vm.count(ss.str().c_str())) {
          species[i]->shadetolerance = vm[ss.str().c_str()].as< double >();
          std::cout << "Using the shade tolerance: " << species[i]->shadetolerance << std::endl;
        }
      }
    }
    if (vm["INTERMEDIATESPECIES"].as< int >() >= 1) {
      std::cout << "Each new tree will be generated with intermediate species parameters (somewhere between White Pine and Eastern Hemlock)" << std::endl;
    }

  }
  catch(std::exception& e)
  {
    std::cout << e.what() << "\n";
  }
  nSeedTotal = initialize_nSeedTotal();
  for (int i=0; i<SMAX; i++) {
    species[i]->nSeed = nSeedTotal / (double)SMAX;
  }
}

