#ifndef SIMOL_OBSERVABLE_HPP
#define SIMOL_OBSERVABLE_HPP

#include <iomanip>
using std::setw;
#include "simol/Tools.hpp"
#include "Statistics.hpp"
#include "simol/system/Particle.hpp"
#include "simol/system/System.hpp"

namespace simol
{
 /// Represents an observable such as the energy or a flux
 /// Statistics such as the mean or the asymptotic variance are estimated on the fly
 /// Manages the ouputs corresponding to this observable
 class Observable
  {
  public:
      // number of steps of the autocorrelation profile used to compute the asymptotic variance
      int decorrelationNbOfSteps_;
      // timestep of the integrator of the dynamics
      double timeStep_;
      // number of steps between two outputs of this observable
      int printPeriodNbOfSteps_;
      // number of nodes in the final plot of the autocorrelation profile (too many nodes can require a very large memory, especially for the chains)
      int nbOfAutocoPts_;
      // Contains the value which is being computed, set to zero at each iteration
      double currentValue_;
      // object used to do statistics on the processed values of the observable
      AutocorrelationStats autocoStats_;
      // memory path to the output file corresponding to this observable 
      string outPath_;
      // output stream for this observable
      ofstream outFlux_;      
      // output stream for the autocorrelation profile of this observable
      shared_ptr<ofstream> outFluxCorrelation_;
    public:
      Observable(Input const& input, int idObs, int decorrelationNbOfSteps0, int nbOfAutocoPts0);
      virtual ~Observable();
      
      const double& timeStep() const;
      const int& decorrelationNbOfSteps() const;
      double decorrelationTime() const;      
      int const& nbOfAutocoPts() const;
      double autocoPtsPeriod() const;
      double printPeriodTime() const;
      const int& printPeriodNbOfSteps() const;
      double time() const;
      
      virtual int nbOfFourier() const
      {return 0;}
      
      virtual int nbOfHermite() const
      {return 0;}
      
      string const& outPath() const;
      ofstream& outFlux();
      ofstream& outFluxCorrelation();
      
      bool doComputeCorrelations() const;

      virtual const double& currentValue() const;
      virtual double& currentValue();
      virtual double lastValue() const;
      virtual double mean() const;
      virtual double asymptoticVariance() const;
      //virtual double stdDev() const;
      virtual double asyvarOfAsyvar() const;
      //virtual double stdDevOfVar() const;

      virtual double correlationAtSpan(int iOfSpan) const;
      virtual double centeredCorrelationAtSpan(int iOfSpan) const;
      virtual double varCorrelationAtSpan(int iOfSpan) const;
      virtual double stdDevCorrelationAtSpan(int iOfSpan) const;

      virtual void append(double value, long int iOfStep);
      virtual void appendCurrent(long int iOfStep);
  
      virtual void updateHamiltonian(System const&){};
      virtual void updateOverdamped(System const&, double /*beta*/){};
      virtual void updateLangevin(System const&, double /*beta*/, double /*gamma*/){};
      virtual void updateBoundaryLangevin(System const&, double /*betaLeft*/, double /*betaRight*/, double /*gamma*/){};

      virtual bool doOutput(long int iOfStep) const;
      //bool doFinal() const;
      virtual void display(long int iOfStep);
      virtual void displayFinalValues(ofstream& out);

      virtual void displayCorrelations(long int iOfStep);
  };
}

#endif