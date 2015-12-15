
#include <cstdlib>
#include <iostream>

#include "multiSystem.hpp"
#include "particleSystem.hpp"
#include "dynamics.hpp"
#include "Vector.hpp"
#include "input.hpp"
#include "galerkin.hpp"

#include <cmath>

#include<yaml-cpp/yaml.h>
#include <getopt.h>

#include "io/CommandLine.hpp"

using std::cout; 
using std::endl; 

int main(int argc, char* argv[])
{
  

  //=====================
  // COMMAND LINE PARSING
  //=====================

  simol::CommandLine cmd(argc,argv);

  //===================
  // INPUT FILE LOADING
  //===================
    
  cout << "Input read in " << cmd.inputFileName() << "...";
  simol::Input input(cmd);

  //============
  // COMPUTATION
  //============
  
  simol::MultiSystem replica(input);
	replica.launch(input);
  
	//simol::Galerkin galerkinSolver(input);
	//galerkinSolver.solve();

  

  std::cout << "Fin de la simulation" << std::endl;

  return EXIT_SUCCESS;
}
