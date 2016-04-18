
#include <cstdlib>
#include <iostream>

#include "simol/statphys/system/System.hpp"
#include "simol/statphys/system/System.hpp"
#include "simol/statphys/dynamics/Dynamics.hpp"
#include "simol/core/linalg/Vector.hpp"
#include "simol/statphys/input/Input.hpp"
#include "simol/statphys/controlVariate/Galerkin.hpp"

#include <cmath>

#include<yaml-cpp/yaml.h>
#include <getopt.h>

#include "simol/core/io/CommandLine.hpp"
#include <time.h>

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
  double totalTime = clock();

  //=====================
  // COMMAND LINE PARSING
  //=====================

  simol::CommandLine cmd(argc, argv);

  //===================
  // INPUT FILE LOADING
  //===================

  cout << "Input read in " << cmd.inputFileName() << endl;
  simol::Input input(cmd);

  //============
  // COMPUTATION
  //============

  if(!input.isGalerkin())
    throw std::invalid_argument("The input file should correspond to a Galerkin computation !");


  simol::BoundaryLangevinGalerkin galerkinSolver(input);
  galerkinSolver.compute();

  displayTime(clock() - totalTime);



  return EXIT_SUCCESS;
}
