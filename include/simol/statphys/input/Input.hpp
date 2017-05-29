#ifndef SIMOL_INPUT_HPP
#define SIMOL_INPUT_HPP

// Disables certain warnings when including the following files
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include<yaml-cpp/yaml.h>
#pragma GCC diagnostic pop

#include <getopt.h>
//#include "simol/core/linalg/Vector.hpp"
#include "simol/core/io/CommandLine.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include "simol/statphys/Tools.hpp"
#include "simol/core/random/RNG.hpp"

namespace simol
{
  double readItem(ifstream& in);
  bool sameLetters(const string& x, const string& y);
  bool isYes(const string& s);
  DVec scanTensor(const string path, vector<int>& dimensions);
  map<double, double> scanMap(const string path);
  string doubleToString(double x);

  class Input
  {
    YAML::Node data;
    string inputPath_;
    ifstream inputFlux_;
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
    //Sinus
    double amplitude() const;
    //DoubleWell
    double height() const;
    double interWell() const;
    //Harmonic
    double potentialCenter() const;
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
    int nbOfAutocoPts() const;
    int nbOfShortAutocoPts() const;
    
    // -- observables --
    string nameOfObs(int idObs) const;
    bool doObservable(int idObs) const;
    bool doKineticEnergy() const;
    bool doPotentialEnergy() const;
    bool doPressure() const;
    bool doDPDE() const;
    bool doInternalEnergy() const;
    bool doInternalTemperature() const;
    bool doVelocity() const;
    bool doForce() const;
    bool doLength() const;
    bool doMidFlow() const;
    bool doSumFlow() const;
    bool doModiFlow() const;
    bool doLagrangeMultiplier() const;
    
    bool doOutThermo() const;
    bool doOutParticles() const;
    bool doXMakeMol() const;
    bool doOutBackUp() const;
    bool doOutChain() const;
    
    bool doFinalLength() const;    
    bool doFinalVelocity() const;
    bool doFinalFlow() const;
    bool doFinalLagrangeMultiplier() const;
    
    bool doOutVelocitiesGenerator() const;
    
    //-- Controle Variate --
    string CVObservable() const;
    bool doCVObservable(int idObs) const;
    string controlVariateName() const;
    string controlVariateCoeffsPath() const;
    bool doGalerkinCV() const;
    bool isGalerkin() const;
    string galerkinElts() const;
    int nbOfQModes() const;
    int nbOfPModes() const;
    bool doGalerkinNonequilibrium() const;
    bool doComputeRef() const;
    bool fitModiFlow() const;
    double omegaHermite() const;
    double integrationStep() const;
    double integrationQMin() const;
    double integrationLength() const;
  };


}

#endif
