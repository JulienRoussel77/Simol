
#include "simol/statphys/dynamics/Hamiltonian.hpp"

namespace simol
{
  //! Constructs a Hamiltonian dynamics (constant energy)
  Hamiltonian::Hamiltonian(Input const& input)
    : Dynamics(input)
  {}
  
  void Hamiltonian::computeThermo(Output& output) const
  {    
    output.temperature() = 2 * output.kineticEnergy() / (output.dimension() * output.nbOfParticles());
    output.totalEnergy() = output.kineticEnergy() + output.potentialEnergy();
    output.pressure() = (2 * output.kineticEnergy() + output.totalVirial()) / (output.dimension() * output.nbOfParticles() * pow(output.latticeParameter(), output.dimension()));
  }
  
  ///
  ///Applies the generator of this dynamics to the basis functions of the CV
  void Hamiltonian::computeGeneratorOnBasis(CVBasis& cvBasis, System const& syst) const
  {
    cvBasis.generatorOnBasisValues_ = DVec::Zero(cvBasis.totalNbOfElts());
    //DVec result = DVec::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < cvBasis.totalNbOfElts(); iOfFunction++)
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        cvBasis.generatorOnBasisValues_(iOfFunction) += dot( syst(iOfParticle).momentum() , cvBasis.basis_->gradientQ(syst, iOfParticle, iOfFunction))
                               + dot( syst(iOfParticle).force() , cvBasis.basis_->gradientP(syst, iOfParticle, iOfFunction));
  }
}
