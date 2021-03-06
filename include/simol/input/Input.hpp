#ifndef SIMOL_INPUT_HPP
#define SIMOL_INPUT_HPP

// Disables certain warnings when including the following files
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include<yaml-cpp/yaml.h>
#pragma GCC diagnostic pop

#include <getopt.h>
#include "simol/CommandLine.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include "simol/Tools.hpp"
#include "simol/RNG.hpp"

namespace simol
{
  double readItem(ifstream& in);
  bool sameLetters(const string& x, const string& y);
  bool isYes(const string& s);
  //DVec scanTensor(const string path, vector<int>& dimensions);
  //map<double, double> scanMap(const string path);
  string doubleToString(double x);

  class Input
  {
    // YAML object managing the fields of the input file
    YAML::Node data;
    // memory address of the input file
    string inputPath_;
    // input flux from the input file
    ifstream inputFlux_;
    // in the case when the system is initialized from a given configuration, this configuration is stored in these vectors
    vector<DVec> initialPositions_, initialMomenta_;
    vector<double> initialInternalEnergies_;
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
    bool doCellMethod() const;
    double mass() const;
    bool doSetting() const;
    bool doRestart() const;
    string restartFileName() const;
    DVec initialPosition(int const& iOfParticle = 0) const;
    DVec initialMomentum(int const& iOfParticle = 0) const;
    double initialInternalEnergy(int const& iOfParticle = 0) const;
    string systemSubtype() const;
    //Chain
    bool isOfFixedVolum() const;
    // NBody
    bool restart() const;

    //-- Time --
    double timeStep() const;
    long int nbOfSteps() const;
    long int thermalizationNbOfSteps() const;
    long int burninNbOfSteps() const;

    //-- Potential --
    string potentialName() const;
    string externalPotentialName() const;
    string pairPotentialName() const;
    string firstPotentialName() const;
    string secondPotentialName() const;
    string galerkinPotentialName() const;
    double interactionRatio() const;
    bool flag() const;
    //Sinus
    double amplitude() const;
    //SpaceSine
    double coupling() const;
    //DoubleWell
    double height() const;
    double interWell() const;
    //Harmonic
    double potentialCenter() const;
    double potentialStiffness() const;
    double potentialAlpha() const;
    double potentialBeta() const;
    double potentialDimension() const;
    // Lennard Jones
    double potentialEpsilon() const;
    double potentialSigma() const;
    double cutOffRatio() const;
    double cutOffRadius() const;
    double splineRatio() const;

    //-- Dynamics --
    string dynamicsName() const;
    bool isConstrained() const;
    //bool isMollified() const;
    double gamma() const;
    double temperature() const;
    bool doMetropolis() const;
    bool doProjectionDPDE() const;
    bool doMTS() const;
    //double temperatureLeft() const;
    //double temperatureRight() const;
    double beta() const;
    //double betaLeft() const;
    //double betaRight() const;
    double deltaTemperature() const;
    double nonEqAmplitude() const;
    double drift() const;
    double tauBending() const;
    double xi() const;
    int seed() const;
    double eta() const;
    double bulkDriving() const;
    double nu() const;
    double flux() const;
    double heatCapacity() const;          // used for DPDE dynamics
    double heatCapacityEinstein() const;  // used for DPDE dynamics
    double kappa() const;         // used for DPDE dynamics
    double einsteinTemperature() const; // used for DPDE dynamics
    int MTSfrequency() const;   // for multiple timestepping, currently DPDE
    double initialInternalTemperature() const; // initial temp. for internal d.o.f. (for equilibration)

    //-- Output --
    string simuTypeName() const;
    //string parametersName() const;
    string outputFolderName() const;
    int decorrelationNbOfSteps() const;
    double decorrelationTime() const;
    int shortDecorrelationNbOfSteps() const;
    double shortDecorrelationTime() const;
    int printPeriodNbOfSteps() const;
    double printPeriodTime() const;
    int printLongPeriodNbOfSteps() const;
    double printLongPeriodTime() const;
    int autocoPtsBinSize() const;
    int nbOfAutocoPts() const;
    int shortAutocoPtsBinSize() const;
    int nbOfShortAutocoPts() const;
    
    // -- observables --
    string nameOfObs(int idObs) const;
    bool doObservable(int idObs) const;
    bool doKineticEnergy() const;
    bool doPotentialEnergy() const;
    bool doTotalEnergy() const;
    bool doPressure() const;
    bool doDPDE() const;
    bool doInternalEnergy() const;
    bool doInternalTemperature() const;
    bool doVelocity() const;
    bool doForce() const;
    bool doLength() const;
    bool doMidFlux() const;
    bool doSumFlux() const;
    bool doModiFlux() const;
    bool doLagrangeMultiplier() const;
    
    bool doOutThermo() const;
    bool doOutParticles() const;
    bool doXMakeMol() const;
    bool doOutBackUp() const;
    bool doOutChain() const;
    
    bool doFinalLength() const;    
    bool doFinalVelocity() const;
    bool doFinalFlux() const;
    bool doFinalLagrangeMultiplier() const;
    bool doFinalChainLagrangeMultiplier() const;
    
  };


}

#endif
