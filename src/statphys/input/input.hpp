#ifndef SIMOL_INPUT_HPP
#define SIMOL_INPUT_HPP

#include<yaml-cpp/yaml.h>
#include <getopt.h>
#include "core/linalg/Vector.hpp"

#include "core/io/CommandLine.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include "tools.hpp"
#include "core/random/RNG.hpp"

namespace simol{

  class Input{
    YAML::Node data;
    string inputPath_;
    ifstream inputFlux_;
    double positionMin_, positionMax_;
    ifstream inputSettings_;
    vector<Vector<double>> initialPositions_, initialMomenta_;
  public:
    Input(CommandLine cmd);
    const string& inputPath() const;
    const std::ifstream& inputFlux() const;
    
    string simuTypeName() const;
    string parametersName() const;
    string outputFolderName() const;
    //const std::shared_ptr<RNG>& rng() const;
    //std::shared_ptr<RNG>& rng();
    
    // Geometry
    int dimension() const;
    double length() const;
    
    // Mesh/Time
    double timeStep() const;
    size_t nbOfIterations() const;
    size_t nbOfThermalIterations() const;
    size_t nbOfBurningIterations() const;
    
    // Physics/System
    string systemName() const;
    size_t nbOfParticles() const;
    double mass() const;
    bool doFileSettings() const;
    string settingsPath() const;
    Vector<double> initialPosition(int const& i=0) const;
    Vector<double> initialMomentum(int const& i=0) const;
    // Physics/Potential
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

    // Physics/Model
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

    // Output
    size_t decorrelationNbOfIterations() const;
    double decorrelationTime() const;
    size_t outputPeriodNbOfIterations() const;
    double outputPeriodTime() const;
    size_t outputProfilePeriodNbOfIterations() const;
    double outputProfilePeriodTime() const;
    int nbOfAutocoPts() const;
    bool doFinalFlow() const;
    bool doFinalVelocity() const;
    
    // Controle Variate
    string controlVariateName() const;
    string controlVariateCoeffsPath() const;
    
    //Galerkin
    bool doGalerkinCV() const;
    bool isGalerkin() const;
    string galerkinElts() const;
    size_t nbOfFourier() const;
    size_t nbOfHermite() const;
    
  };


}

#endif
