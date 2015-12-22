
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
#include <time.h>

using std::cout; 
using std::endl; 

int main(int argc, char* argv[])
{
  double totalTime = clock();

  //=====================
  // COMMAND LINE PARSING
  //=====================

  simol::CommandLine cmd(argc,argv);

  //===================
  // INPUT FILE LOADING
  //===================
    
  cout << "Input read in " << cmd.inputFileName() << endl;
  simol::Input input(cmd);

  //============
  // COMPUTATION
  //============
  
  simol::MultiSystem replica(input);
	replica.launch(input);
  
	//simol::Galerkin galerkinSolver(input);
	//galerkinSolver.compute();

  

  std::cout << "Fin de la simulation" << std::endl;
	int nbSeconds = (clock() - totalTime) / CLOCKS_PER_SEC;
	int nbMinutes = nbSeconds / 60;
	int nbHours = nbMinutes / 60;
	int nbDays = nbHours / 24;
	cout << "Temps de calcul total : " << nbDays << " j " <<  nbHours%24  << " h " << nbMinutes%60 << " m " << nbSeconds%60 << " s " << endl;

  return EXIT_SUCCESS;
}
