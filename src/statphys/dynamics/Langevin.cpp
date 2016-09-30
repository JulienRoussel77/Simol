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
  ///Print de name of the dynamics
  void Langevin::printName() const
  {
    std::cout << "Dynamics = Langevin" << std::endl;
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
  void Langevin::computeGeneratorOnBasis(CVBasis& cvBasis, vector<Particle*> const& configuration) const
  {
    int nbOfParticles = (int)configuration.size();
    //cout << "generatorOn(const Langevin& dyna, S const& syst, const ControlVariate& controlVariate)" << endl;
    //Vector<double> result = Vector<double>::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < cvBasis.nbOfFunctions(); iOfFunction++)
      for (int iOfParticle = 0; iOfParticle < nbOfParticles; iOfParticle++)
      {
        cvBasis.generatorOnBasisValues_(iOfFunction) += dot(configuration[iOfParticle]->momentum() , cvBasis.basis_->gradientQ(configuration, iOfParticle, iOfFunction))
                               + dot(configuration[iOfParticle]->force() , cvBasis.basis_->gradientP(configuration, iOfParticle, iOfFunction))
                               + gamma() * (- dot(configuration[iOfParticle]->momentum() , cvBasis.basis_->gradientP(configuration, iOfParticle, iOfFunction))
                                          + cvBasis.basis_->laplacianP(configuration, iOfParticle, iOfFunction) / beta() );
      }
  }

}
