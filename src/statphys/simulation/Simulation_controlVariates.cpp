#include "Simulation.hpp"

namespace simol {
  
  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(Dynamics const& /*dyna*/, System const& syst, Output& /*output*/, int /*iOfIteration*/)
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
    
  void updateAllControlVariates(const Hamiltonian& /*dyna*/, System const& syst, Output& /*output*/, int /*iOfIteration*/)
  {}
  
  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(const Langevin& dyna, System const& syst, Output& output, int iOfIteration)
  {
    //cout << "updateAllControlVariates(const Langevin& dyna, S const& syst, Output& output, int iOfIteration)" << endl;
    Vector<double> q = syst.getParticle(0).position();
    Vector<double> p = syst.getParticle(0).momentum();
    Vector<double> generatorOnBasis;
    //generatorOnBasis = generatorOn(dyna, syst, output.velocityCV());
    generatorOnBasis = output.velocityCV().generatorLangevin(syst.configuration(), dyna.beta(), dyna.gamma());
    output.velocityCV().update(p(0), generatorOnBasis, syst.configuration(), iOfIteration);
    if (output.doOutput(iOfIteration))
      output.displayGeneratorOnBasis(output.outVelocitiesGenerator(), syst.configuration(), output.velocityCV(), iOfIteration*dyna.timeStep());
  }
  
  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(const BoundaryLangevin& dyna, System const& syst, Output& output, int iOfIteration)
  {
    Vector<double> generatorOnBasis;
    
    //generatorOnBasis = generatorOn(dyna, syst, output.midFlowCV());
    generatorOnBasis = output.midFlowCV().generatorBoundarylangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma());
    output.midFlowCV().update(output.energyMidFlow(), generatorOnBasis, syst.configuration(), iOfIteration);
    //generatorOnBasis = generatorOn(dyna, syst, output.sumFlowCV());
    generatorOnBasis = output.sumFlowCV().generatorBoundarylangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma());
    output.sumFlowCV().update(output.energySumFlow(), generatorOnBasis, syst.configuration(), iOfIteration);
  }

}
