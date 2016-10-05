#include "simol/statphys/dynamics/Langevin.hpp"

namespace simol
{
  ///Constructor for the Langevin Dynamics with thermostats everywhere
  Langevin::Langevin(Input const& input):
    LangevinBase(input)
  {
    galerkin_ = createLangevinGalerkin(input);
  }


  ///
  ///Returns the amplitude of the brownian motion
  double Langevin::sigma() const
  {
    return sqrt(2 * gamma_ / beta_);
  }

  /*///After refers to the fact that this step comes after updating the forces
  ///Proceeds to the second half of the Verlet scheme, then integrate the O-U analytically
  void Langevin::updateAfter(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1 - pow(alpha, 2)) / beta_ * particle.mass()) * rng_->gaussian();
  }*/
  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  void Langevin::computeGeneratorOnBasis(CVBasis& cvBasis, System const& syst) const
  {
    //cout << "generatorOn(const Langevin& dyna, S const& syst, const ControlVariate& controlVariate)" << endl;
    //Vector<double> result = Vector<double>::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < cvBasis.nbOfFunctions(); iOfFunction++)
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      {
        cvBasis.generatorOnBasisValues_(iOfFunction) += dot(syst(iOfParticle).momentum() , cvBasis.basis_->gradientQ(syst, iOfParticle, iOfFunction))
                               + dot(syst(iOfParticle).force() , cvBasis.basis_->gradientP(syst, iOfParticle, iOfFunction))
                               + gamma() * (- dot(syst(iOfParticle).momentum() , cvBasis.basis_->gradientP(syst, iOfParticle, iOfFunction))
                                          + cvBasis.basis_->laplacianP(syst, iOfParticle, iOfFunction) / beta() );
      }
  }

}
