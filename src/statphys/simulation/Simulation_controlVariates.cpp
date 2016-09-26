#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{

  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(Dynamics const& /*dyna*/, System const& /*syst*/, Output& /*output*/, long int /*iOfStep*/)
  {throw std::invalid_argument("updateAllControlVariates: Function undefined");}

  void updateAllControlVariates(const Hamiltonian& /*dyna*/, System const& /*syst*/, Output& /*output*/, long int /*iOfStep*/)
  {}
  
  void updateAllControlVariates(const Overdamped& /*dyna*/, System const& /*syst*/, Output& /*output*/, long int /*iOfStep*/)
  {}

  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(const Langevin& dyna, System const& syst, Output& output, long int iOfStep)
  {
    //cout << "updateAllControlVariates(const Langevin& dyna, S const& syst, Output& output, long int iOfStep)" << endl;
    Vector<double> q = syst.getParticle(0).position();
    Vector<double> p = syst.getParticle(0).momentum();
    
    output.obsVelocity().updateLangevin(syst.configuration(), dyna.beta(), dyna.gamma());
    output.obsVelocity().append(p(0), iOfStep);
    
    /*if (!output.obsVelocity().hasControlVariate())    
      output.obsVelocity().update(p(0), iOfStep);
    else{
      output.obsVelocity().computeGeneratorLangevin(syst.configuration(), dyna.beta(), dyna.gamma());
      //Vector<double> generatorOnBasis;
      //generatorOnBasis = output.obsVelocity().generatorLangevin(syst.configuration(), dyna.beta(), dyna.gamma());
      output.obsVelocity().update(p(0), syst.configuration(), iOfStep);
      if (output.doOutput(iOfStep))
        output.displayGeneratorOnBasis(output.outVelocitiesGenerator(), syst.configuration(), output.obsVelocity(), iOfStep * dyna.timeStep());
    }*/
  }

  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(const BoundaryLangevin& dyna, System const& syst, Output& output, long int iOfStep)
  {
    output.obsMidFlow().updateBoundaryLangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma());
    //output.obsMidFlow().append(output.energyMidFlow(), iOfStep);
    output.obsMidFlow().appendCurrent(iOfStep);
    
    output.obsSumFlow().updateBoundaryLangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma());
    //output.obsSumFlow().append(output.energySumFlow(), iOfStep);
    output.obsSumFlow().appendCurrent(iOfStep);
    
    output.obsModiFlow().updateBoundaryLangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma());
    //output.obsModiFlow().append(output.energyModiFlow(), iOfStep);
    output.obsModiFlow().appendCurrent(iOfStep);

    /*if (!output.obsMidFlow().hasControlVariate())
      output.obsMidFlow().append(output.energyMidFlow(), iOfStep);
    else
    {
      output.obsMidFlow().computeGeneratorBoundarylangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma());
      output.obsMidFlow().append(output.energyMidFlow(), syst.configuration(), iOfStep);
    }
    
    if (!output.obsSumFlow().hasControlVariate())
      output.obsSumFlow().append(output.energySumFlow(), iOfStep);
    else
    {
      output.obsSumFlow().computeGeneratorBoundarylangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma());
      output.obsSumFlow().append(output.energySumFlow(), syst.configuration(), iOfStep);
    }
    
    if (!output.obsModiFlow().hasControlVariate())
      output.obsModiFlow().append(output.energyModiFlow(), iOfStep);
    else
    {
      output.obsModiFlow().computeGeneratorBoundarylangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma());
      output.obsModiFlow().append(output.energyModiFlow(), syst.configuration(), iOfStep);
    }*/
    //*(output.outTest_) << iOfStep << " " << output.obsModiFlowCV().statsObservable_.statsValues_.lastValue() * output.obsModiFlowCV().statsObservable_.statsRefValues_.lastValue() << endl;
  }

}
