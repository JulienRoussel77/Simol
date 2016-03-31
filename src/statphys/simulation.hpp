#ifndef SIMOL_SIMULATION_HPP
#define SIMOL_SIMULATION_HPP

#include "tools.hpp"
#include "system.hpp"
#include "chain/chain.hpp"
#include "controlVariate.hpp"

#include "dynamics/BoundaryLangevin.hpp"
#include "dynamics/Hamiltonian.hpp"
#include "dynamics/Overdamped.hpp"
#include "dynamics/Langevin.hpp"

namespace simol {

  template<typename D, typename S>
  class Simulation
  {
    int dimension_;
    std::shared_ptr<RNG> rng_;
    S system_;
    D dynamics_;
    Output output_;
    //ControlVariate* controlVariate_;
  public:
    Simulation(Input& input);
    //virtual ~Simulation();

    void launch();
  };

    template<class D, class S>
  Simulation<D,S>::Simulation(Input& input):
    dimension_(input.dimension()),
    rng_(std::make_shared<RNG>(RNG(input.seed(), input.dimension()))),
    system_(input),
    dynamics_(input),
    output_(input)
 {
   //cout << "Simulation::Simulation(Input& input)" << endl;
   output_.verbose() = 1;

   system_.rng() = rng_;
   dynamics_.rng() = rng_;

   output_.setControlVariates(input, system_.potential(), dynamics_.galerkin());
 }

  template<class D, class S>
  void Simulation<D,S>::launch()
  {
    //cout << "Simulation::launch" << endl;
    launchSimu(dynamics_, system_, output_);
  }



  // --------------- Declaration and implementation of external functions -------------------



  template<class D, class S>
  void initializeMomenta(const D& dyna, S& syst);
  template<class D, class S>
  void initializeSystem(D& dyna, S& syst);
  template<class D, class S>
  void simulate(D& dyna, S& syst);
  template<class D, class S>
  void updateAllControlVariates(const D& dyna, const S& syst, Output& output, size_t iOfIteration);
  template<class D, class S>
  Vector<double> generatorOn(const D& dyna, const S& syst, const ControlVariate& controlVariate);
  template <class D, class S>
  void computeOutput(D const& dyna, S const& syst, Output& output, size_t iOfIteration);
  template <class S>
  void writeOutput(S const& syst, Output& output, size_t iOfIteration);
  template <class D, class S>
  void writeFinalOutput(D const& dyna, S const& syst, Output& output);
  template <class D, class S>
  void updateAllControlVariates(D const& dyna, S const& syst, Output& output, size_t iOfIteration);
  template <class D, class S>
  void simulate(D& dyna, S& syst);

  //Isolated
  template <class D>
  void initializeSystem(Dynamics& dyna, Isolated& syst);
  template <class D>
  void computeOutput(D const& dyna, const Isolated& syst, Output& output, size_t iOfIteration);
  template <class S>
  void writeOutput(S const& syst, Output& output, size_t iOfIteration);
  template <class D, class S>
  void writeFinalOutput(D const& dyna, S const& syst, Output& output);

  //Fluid
  template <class D, class S>
  void writeFinalOutput(D const& dyna, S const& syst, Output& output);

  //Hamiltonian
  template <class S>
  Vector<double> generatorOn(const Hamiltonian& dyna, S const& syst, ControlVariate const& controlVariate);
  template <class S>
  void updateAllControlVariates(const Hamiltonian& dyna, S const& syst, Output& output, size_t iOfIteration);

  //Langevin
  template <class S>
  Vector<double> generatorOn(const Langevin& dyna, S const& syst, const ControlVariate& controlVariate);
  template <class S>
  void updateAllControlVariates(const Langevin& dyna, S const& syst, Output& output, size_t iOfIteration);

  //Overdamped
  template <class S>
  Vector<double> generatorOn(const Overdamped& dyna, const System& syst, const ControlVariate& controlVariate);

  //Chains
  template <>
  void initializeMomenta(const BoundaryLangevin& dyna, Chain& syst);
  template <>
  void initializeSystem(BoundaryLangevin& dyna, BiChain& syst);
  template <>
  void initializeSystem(BoundaryLangevin& dyna, TriChain& syst);
  template <>
  void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output);
  template <>
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output);
  template <>
  void simulate(BoundaryLangevin& dyna, Chain& syst);
  template <class S>
  Vector<double> generatorOn(const BoundaryLangevin& dyna, S const& syst, const ControlVariate& controlVariate);

  //DPDE
  template <class S>
  void simulate(DPDE& dyna, S& syst);
  template<>
  void initializeSystem(DPDE& dyna, Isolated& syst);
  template<>
  void computeOutput(const DPDE& dyna, const Isolated& syst, Output& output, size_t iOfIteration);
  template <class S>
  void writeOutput(DPDE const& dyna, S const& syst, Output& output, size_t iOfIteration);

  //General
  template<class D, class S>
  void launchSimu(D& dyna, S& syst, Output& output);











  // -----------------------------Implementation des templates----------------------------------



  ///
  ///Initializes all the momenta, assuming an affine temperature profile
  template <class D, class S>
  void initializeMomenta(D const& /*dyna*/, S& /*syst*/)
  {throw std::invalid_argument("initializeMomenta : Function undefined");}

  template <class D, class S>
  void initializeSystem(D& /*dyna*/, S& /*syst*/)
  {throw std::invalid_argument("initializeSystem : Function undefined");}

   ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  template <class D, class S>
  void updateAllControlVariates(D const& dyna, S const& syst, Output& output, size_t iOfIteration)
  {throw std::invalid_argument("updateAllControlVariates: Function undefined");}

    /*Vector<double> q = syst.getParticle(0).position();
    Vector<double> p = syst.getParticle(0).momentum();
    Vector<double> qEnd = syst.getParticle(syst.nbOfParticles()-1).position();
    Vector<double> generatorOnBasis;
    generatorOnBasis = generatorOn(dyna, syst, output.velocityCV());
    output.velocityCV().update(p(0), generatorOnBasis, syst.configuration(), iOfIteration);
    generatorOnBasis = generatorOn(dyna, syst, output.forceCV());
    output.forceCV().update(syst.force(q)(0), generatorOnBasis, syst.configuration(), iOfIteration);
    generatorOnBasis = generatorOn(dyna, syst, output.lengthCV());
    output.lengthCV().update(qEnd(0), generatorOnBasis, syst.configuration(), iOfIteration);
    generatorOnBasis = generatorOn(dyna, syst, output.midFlowCV());
    output.midFlowCV().update(output.energyMidFlow(), generatorOnBasis, syst.configuration(), iOfIteration);
    generatorOnBasis = generatorOn(dyna, syst, output.sumFlowCV());
    output.sumFlowCV().update(output.energySumFlow(), generatorOnBasis, syst.configuration(), iOfIteration);
  }*/

  ///
  ///Applies the generator of this dynamics to the basis functions of the CV
  template <class D, class S>
  Vector<double> generatorOn(D const& /*dyna*/, S const& /*syst*/, const ControlVariate& /*controlVariate*/)
  {
    throw std::invalid_argument("GeneratorOn not defined in the general case");
  }

  template <class D, class S>
  void computeOutput(D const& dyna, S const& syst, Output& output, size_t iOfIteration)
  {
    if (output.verbose() > 0)
    {
      output.kineticEnergy() = 0;
      output.potentialEnergy() = 0;
      //Calcul de la température et de l'énergie
      for (const auto& particle : syst.configuration())
      {
        //Particle& particle = syst.getParticle(iOfParticle);
        output.kineticEnergy() += particle.kineticEnergy();
        output.potentialEnergy() += particle.potentialEnergy();
      }
    }
    // In the case of the trichain we add the potential of the wall interaction
    output.potentialEnergy() += syst.boundaryPotEnergy();
    syst.computeProfile(output, dyna, iOfIteration);
    updateAllControlVariates(dyna, syst, output, iOfIteration);
  }

  // TO DO : cette fonction ne sert pas a grand chose... ce test peut (DOIT) etre fait dans output.cpp !!
  template <class D, class S>
  void writeOutput(D const& /*dyna*/, S const& syst, Output& output, size_t iOfIteration)
  {
    if (output.verbose() > 0 && output.doOutput(iOfIteration))// && iOfIteration >= 100)
      output.display(syst.configuration(), iOfIteration);
  }

  template <class D, class S>
  void writeFinalOutput(D const& dyna, S const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
  }

  template <class D, class S>
  void simulate(D& dyna, S& syst)
  {
    for (auto&& particle : syst.configuration())
      dyna.updateBefore(particle);

    syst.computeAllForces(dyna);

    for (auto&& particle : syst.configuration())
      dyna.updateAfter(particle);
  }





 //##################### ISOLATED ######################

  // TODO : beta shouldn't be called here
  template <class D>
  void initializeSystem(D& dyna, Isolated& syst)
  {
    syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
    syst.getParticle(0).position(0) = syst.drawPotLaw(dyna.beta());
  }

  template <class D>
  void computeOutput(D const& dyna, const Isolated& syst, Output& output, size_t iOfIteration)
  {
    //cout << "computeOutput(D const& dyna, const Isolated& syst, Output& output, size_t iOfIteration)" << endl;
    //Calcul de la température et de l'énergie
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    updateAllControlVariates(dyna, syst, output, iOfIteration);
  }

  template <class D>
  void writeFinalOutput(D const& dyna, Isolated const& syst, Output& output)
  {
    //double time = dyna.timeStep() * dyna.nbOfIterations();

    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
    output.displayFinalVelocity(dyna.temperature(), dyna.externalForce(0), output.velocityCV_->nbOfFourier(), output.velocityCV_->nbOfHermite());
  }


  //################## FLUID ########################

  template <class D>
  void writeFinalOutput(D const& dyna, Fluid const& syst, Output& output)
  {
    //double time = dyna.timeStep() * dyna.nbOfIterations();

    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
    //output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature());
  }

  //##################### HAMILTONIAN ####################"

  template <>
  void initializeSystem(Hamiltonian& /*dyna*/, Isolated& /*syst*/)
  {
    //syst.getParticle(0).momentum().fill(0);
    //syst.getParticle(0).position().fill(0);
  }

  ///
  ///Applies the generator of this dynamics to the basis functions of the CV
  template <class S>
  Vector<double> generatorOn(const Hamiltonian& dyna, S const& syst, ControlVariate const& controlVariate)
  {
    Vector<double> result = Vector<double>::Zero(controlVariate.nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate.nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        result(iOfFunction) += syst.getParticle(iOfParticle).momentum().dot(controlVariate.gradientQ(syst.configuration(), iOfParticle, iOfFunction))
        + syst.force(syst.getParticle(iOfParticle).position()).dot(controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction))
        + dyna.externalForce().dot(controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction));
    return result;
  }

  template <class S>
  void updateAllControlVariates(const Hamiltonian& /*dyna*/, S const& /*syst*/, Output& /*output*/, size_t /*iOfIteration*/)
  {}

  template <>
  void writeFinalOutput(Hamiltonian const& dyna, Isolated const& syst, Output& output)
  {
    //double time = dyna.timeStep() * dyna.nbOfIterations();

    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
    output.displayFinalVelocity(0, dyna.externalForce(0), output.velocityCV_->nbOfFourier(), output.velocityCV_->nbOfHermite());
  }


  //################# LANGEVIN ###########################

  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  template <class S>
  Vector<double> generatorOn(const Langevin& dyna, S const& syst, const ControlVariate& controlVariate)
  {
    //cout << "generatorOn(const Langevin& dyna, S const& syst, const ControlVariate& controlVariate)" << endl;
    Vector<double> result = Vector<double>::Zero(controlVariate.nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate.nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      {
        result(iOfFunction) += dot(syst.getParticle(iOfParticle).momentum(), controlVariate.gradientQ(syst.configuration(), iOfParticle, iOfFunction))
          + dot(syst.getParticle(iOfParticle).force(), controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction))
          + dyna.gamma() * (- dot(syst.getParticle(iOfParticle).momentum(), controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction))
              + controlVariate.laplacianP(syst.configuration(), iOfParticle, iOfFunction) / dyna.beta() );
      }
    return result;
  }

  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  template <class S>
  void updateAllControlVariates(const Langevin& dyna, S const& syst, Output& output, size_t iOfIteration)
  {
    //cout << "updateAllControlVariates(const Langevin& dyna, S const& syst, Output& output, size_t iOfIteration)" << endl;
    Vector<double> q = syst.getParticle(0).position();
    Vector<double> p = syst.getParticle(0).momentum();
    Vector<double> generatorOnBasis;
    generatorOnBasis = generatorOn(dyna, syst, output.velocityCV());
    output.velocityCV().update(p(0), generatorOnBasis, syst.configuration(), iOfIteration);
    if (output.doOutput(iOfIteration))
      output.displayGeneratorOnBasis(output.outVelocitiesGenerator_, syst.configuration(), output.velocityCV(), iOfIteration*dyna.timeStep());

    /*generatorOnBasis = generatorOn(output.forceCV(), configuration);
    output.forceCV().update(potential_->derivative(q)(0), generatorOnBasis, configuration, iOfIteration);
    generatorOnBasis = generatorOn(output.lengthCV(), configuration);
    output.lengthCV().update(q(0), generatorOnBasis, configuration, iOfIteration);*/
    //cout << "end updateAllControlVariates" << endl;
  }

//################### OVERDAMPED ########################

  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  template<class S>
  Vector<double> generatorOn(const Overdamped& dyna, S const& syst, const ControlVariate& controlVariate)
  {
    Vector<double> result = Vector<double>::Zero(controlVariate.nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate.nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        result(iOfFunction) += controlVariate.laplacianQ(syst.configuration(), iOfParticle, iOfFunction) / dyna.beta()
          + dot(syst.force(syst.getParticle(iOfParticle).position()), controlVariate.gradientQ(syst.configuration(), iOfParticle, iOfFunction));
    return result;
  }

  //################### CHAINS #############################

  ///
  ///Initializes all the momenta, assuming an affine temperature profile
  template <>
  void initializeMomenta(const BoundaryLangevin& dyna, Chain& syst)
  {
    for (size_t iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      Particle& particle = syst.getParticle(iOfParticle);
      double tempi_ = dyna.temperatureLeft() + (iOfParticle + .5) * (dyna.temperatureRight() - dyna.temperatureLeft()) / syst.nbOfParticles();
      particle.momentum() = sqrt(tempi_ / particle.mass()) * syst.rng()->gaussian();
    }
  }

  template <>
  void initializeSystem(BoundaryLangevin& dyna, BiChain& syst)
  {
    cout << "Initialization of the system...";cout.flush();
    double alpha, localTemp, localDist;
    //ofstream outTest("test.txt");
    for (size_t iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      alpha = iOfParticle / (double) syst.nbOfParticles();
      localTemp = (1-alpha) * dyna.temperatureLeft() + alpha * dyna.temperatureRight();
      //localBending = dyna.computeMeanPotLaw(1/localTemp);
      localDist = syst.drawPotLaw(1/localTemp);
      //outTest << iOfParticle << " " << localBending << endl;
      //cout << "bending = " << localBending << " / mean = " << dyna.computeMeanPotLaw(1/localTemp) << endl;

      double prevPosition = (iOfParticle>0)?syst.getParticle(iOfParticle-1).position(0):0;


      syst.getParticle(iOfParticle).position(0) = prevPosition + localDist;
      //cout << -position2 << " + 2 * " << position1 << " + " << localBending << " = " << syst.getParticle(iOfParticle).position(0) << endl;
      syst.getParticle(iOfParticle).momentum() = syst.drawMomentum(1/localTemp, syst.getParticle(iOfParticle).mass());

      dyna.initializeCountdown(syst.getParticle(iOfParticle));
    }

    cout << "Done ! / Thermalization...";cout.flush();

    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfThermalIterations(); ++iOfIteration)
    {
      syst.thermalize(dyna);
    }

    cout << "Done ! / Burning...";cout.flush();

    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfBurningIterations(); ++iOfIteration)
    {
      simulate(dyna, syst);
    }
    cout << "Done !" << endl;
  }

  template <>
  void initializeSystem(BoundaryLangevin& dyna, TriChain& syst)
  {
    cout << "Initialization of the system...";cout.flush();
    double alpha, localTemp, localBending;
    //ofstream outTest("test.txt");
    for (size_t iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      alpha = iOfParticle / (double) syst.nbOfParticles();
      localTemp = (1-alpha) * dyna.temperatureLeft() + alpha * dyna.temperatureRight();
      //localBending = dyna.computeMeanPotLaw(1/localTemp);
      localBending = syst.drawPotLaw(1/localTemp);
      //outTest << iOfParticle << " " << localBending << endl;
      //cout << "bending = " << localBending << " / mean = " << dyna.computeMeanPotLaw(1/localTemp) << endl;

      double position1 = (iOfParticle>0)?syst.getParticle(iOfParticle-1).position(0):0;
      double position2 = (iOfParticle>1)?syst.getParticle(iOfParticle-2).position(0):0;


      syst.getParticle(iOfParticle).position(0) = -position2 + 2 * position1 + localBending;
      //cout << -position2 << " + 2 * " << position1 << " + " << localBending << " = " << syst.getParticle(iOfParticle).position(0) << endl;
      syst.getParticle(iOfParticle).momentum() = syst.drawMomentum(1/localTemp, syst.getParticle(iOfParticle).mass());

      dyna.initializeCountdown(syst.getParticle(iOfParticle));
    }

    cout << "Done ! / Thermalization...";cout.flush();

    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfThermalIterations(); ++iOfIteration)
    {
      syst.thermalize(dyna);
    }

    cout << "Done ! / Burning...";cout.flush();

    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfBurningIterations(); ++iOfIteration)
    {
      simulate(dyna, syst);
    }
    cout << "Done !" << endl;
  }

  template <>
  void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output)
  {
    //double time = dyna.timeStep() * dyna.nbOfIterations();

    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature());
  }

  template <>
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output)
  {
    //double time = dyna.timeStep() * dyna.nbOfIterations();

    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature(), dyna.tauBending(), dyna.xi());
  }


  template <>
  void simulate(BoundaryLangevin& dyna, Chain& syst)
  {
    //for (auto&& particle : configuration_)
      //dyna.updateBefore(particle);
    //for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    //  dyna.updateBefore(getParticle(iOfParticle));

    for (auto&& particle:syst.configuration())
      dyna.updateBefore(particle);

    syst.computeAllForces(dyna);

    for (auto&& particle:syst.configuration())
      dyna.updateAfter(particle);

    //for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    //  dyna.updateAfter(getParticle(iOfParticle));

    dyna.updateOrsteinUhlenbeck(syst.getParticle(0), dyna.betaLeft());
    dyna.updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), dyna.betaRight());

    if (dyna.doMomentaExchange())
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles()-1; iOfParticle++)
        dyna.updateMomentaExchange(syst.getParticle(iOfParticle), syst.getParticle(iOfParticle+1));
  }


  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  template <class S>
  Vector<double> generatorOn(const BoundaryLangevin& dyna, S const& syst, const ControlVariate& controlVariate)
  {
    size_t nbOfParticles = syst.nbOfParticles();
    Vector<double> result = Vector<double>::Zero(controlVariate.nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate.nbOfFunctions(); iOfFunction++)
    {
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      result(iOfFunction) += dot(syst.getParticle(iOfParticle).momentum(), controlVariate.gradientQ(syst.configuration(), iOfParticle, iOfFunction))
        + dot(syst.getParticle(iOfParticle).force(), controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction));
        //if(false)
      result(iOfFunction) += dyna.gamma() * (- dot(syst.getParticle(0).momentum(), controlVariate.gradientP(syst.configuration(), 0, iOfFunction))
      + controlVariate.laplacianP(syst.configuration(), 0, iOfFunction) / dyna.betaLeft()
      - dot(syst.getParticle(nbOfParticles-1).momentum(), controlVariate.gradientP(syst.configuration(), nbOfParticles-1, iOfFunction))
      + controlVariate.laplacianP(syst.configuration(), nbOfParticles-1, iOfFunction) / dyna.betaRight());
    }
    return result;
  }

  //-------------- DPDE -------------

  template <class S>
  void simulate(DPDE& dyna, S& syst)
  {
    //-- Verlet part -- TO DO : a completer...
    for (auto&& particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces(dyna);
    for (auto&& particle : syst.configuration())
      dyna.verletSecondPart(particle);
    //-- fluctuation/dissipation --
    //for (auto&& particle : syst.configuration())
    for (size_t i=0; i < syst.nbOfParticles(); i++)  // version explicite de la ligne cachee ci dessus, utile pour faire des boucles sur les couples !
      dyna.energyReinjection(syst.getParticle(i));  // integration de p avec gamma fixe + reinjection
  }

  template<>
  void initializeSystem(DPDE& dyna, Isolated& syst)
  {
    syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
    syst.getParticle(0).position(0) = 0;
    syst.getParticle(0).internalEnergy() = 1;  // TO DO : il faudra ici tirer selon la bonne loi ?
  }

  template<>
  void computeOutput(const DPDE& /*dyna*/, const Isolated& syst, Output& output, size_t /*iOfIteration*/)
  {
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    output.internalEnergy() = syst.getParticle(0).internalEnergy();
  }

  template <class S>
  void writeOutput(DPDE const& /*dyna*/, S const& syst, Output& output, size_t iOfIteration)
  {
    if (output.verbose() > 0 && output.doOutput(iOfIteration))
      output.display_DPDE(syst.configuration(), iOfIteration);
  }





      //------------------------------ MAIN FUNCTION ------------------------------------

  /*void printClass(System& syst);
  void printClass(Isolated& syst);
  void printClass(Langevin& syst);
  void printClass(Dynamics& syst);*/

  template <class D, class S>
  void launchSimu(D& dyna, S& syst, Output& output)
  {
    cout << "Estimated time : " << 3.5 * syst.nbOfParticles()/1024. * dyna.nbOfIterations() / 1e6 << " hours" << endl;

    //---- initialization (including burn-in) -----
    initializeSystem(dyna, syst);
    syst.computeAllForces(dyna);  // TO DO : a mettre dans la fonction d'initialisation...

    //---- actual iterations -----
    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfIterations(); ++iOfIteration)
      {
        //--- display progress every time 10% of simulation elapsed ---
        if ((10*iOfIteration) % dyna.nbOfIterations() == 0)
          cout << "---- Run " << (100 * iOfIteration) / dyna.nbOfIterations() << " % completed ----" << endl;

        //--- write outputs if required ----
        computeOutput(dyna, syst, output, iOfIteration);
        writeOutput(dyna, syst, output, iOfIteration);

        //---- update the system by the numerical integration ---
        simulate(dyna, syst);
      }

    //--- write final outputs ----
    writeFinalOutput(dyna, syst, output);
  }


}

#endif
