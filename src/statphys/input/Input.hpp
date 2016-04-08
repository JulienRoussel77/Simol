#ifndef SIMOL_INPUT_HPP
#define SIMOL_INPUT_HPP

#include<yaml-cpp/yaml.h>
#include <getopt.h>
#include "core/linalg/Vector.hpp"
#include "core/io/CommandLine.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include "Tools.hpp"
#include "core/random/RNG.hpp"

namespace simol{

  class Input{
    YAML::Node data;
    string inputPath_;
    ifstream inputFlux_;
    ifstream inputSettings_;
    vector<Vector<double>> initialPositions_, initialMomenta_;
  public:
    Input(CommandLine cmd);
    const string& inputPath() const;
    const std::ifstream& inputFlux() const;
    
    //-- System --
    int dimension() const;
    double length() const;
    double latticeParameter() const;
    string systemName() const;
    int nbOfParticles() const;
    int nbOfParticlesPerDimension() const;
    double mass() const;
    bool doFileSettings() const;
    string settingsPath() const;
    Vector<double> initialPosition(int const& i=0) const;
    Vector<double> initialMomentum(int const& i=0) const;
    
    //-- Time --
    double timeStep() const;
    int nbOfIterations() const;
    int nbOfThermalIterations() const;
    int nbOfBurningIterations() const;
    
    //-- Potential --
    string potentialName() const;
    //Sinus
    double amplitude() const;
    //DoubleWell
    double height() const;
    double interWell() const;
    //Harmonic
    double potentialStiffness() const;
    double potentialAlpha() const;
    double potentialBeta() const;
    // Lennard Jones
    double potentialEpsilon() const;
    double potentialSigma() const;
    double cutOffRatio() const;
    double splineRatio() const;
    
    //-- Dynamics --
    string dynamicsName() const;
    double gamma() const;
    double temperature() const;
    double temperatureLeft() const;
    double temperatureRight() const;
    double beta() const;
    double betaLeft() const;
    double betaRight() const;
    double externalForce() const;
    double tauBending() const;
    double xi() const;
    int seed() const;
    double eta() const;
    double heatCapacity() const;  // used for DPDE dynamics

    //-- Output --
    string simuTypeName() const;
    string parametersName() const;
    string outputFolderName() const;
    int decorrelationNbOfIterations() const;
    double decorrelationTime() const;
    int outputPeriodNbOfIterations() const;
    double outputPeriodTime() const;
    int outputProfilePeriodNbOfIterations() const;
    double outputProfilePeriodTime() const;
    int nbOfAutocoPts() const;
    bool doFinalFlow() const;
    bool doFinalVelocity() const;
    
    //-- Controle Variate --
    string controlVariateName() const;
    string controlVariateCoeffsPath() const;
    bool doGalerkinCV() const;
    bool isGalerkin() const;
    string galerkinElts() const;
    int nbOfFourier() const;
    int nbOfHermite() const;
    
  };


}

#endif
