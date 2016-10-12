#include "simol/statphys/dynamics/Langevin.hpp"

namespace simol
{
  ///Constructor for the Langevin Dynamics with thermostats everywhere
  Langevin::Langevin(Input const& input):
    LangevinBase(input)
  {
    //galerkin_ = createLangevinGalerkin(input);
    galerkin_ = createGalerkin(input);
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
    cvBasis.generatorOnBasisValues_ = DVec::Zero(cvBasis.totalNbOfElts());
    for (int iOfFunction = 0; iOfFunction < cvBasis.totalNbOfElts(); iOfFunction++)
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      {
        cvBasis.generatorOnBasisValues_(iOfFunction) += dot(syst(iOfParticle).momentum() , cvBasis.basis_->gradientQ(syst, iOfParticle, iOfFunction))
                               + dot(syst(iOfParticle).force() , cvBasis.basis_->gradientP(syst, iOfParticle, iOfFunction))
                               + gamma() * (- dot(syst(iOfParticle).momentum() , cvBasis.basis_->gradientP(syst, iOfParticle, iOfFunction))
                                          + cvBasis.basis_->laplacianP(syst, iOfParticle, iOfFunction) / beta() );
      }
      
    //ofstream tempOut("aaaaaaaaaaaaaaaaaaaa.txt", std::ofstream::app);
    //tempOut << syst(0).position(0) << " " << syst(0).momentum(0) << " " << syst(0).force(0) << " " << dot(*cvBasis.cvCoeffs_, cvBasis.generatorOnBasisValues_) << endl;
  }
  


}
