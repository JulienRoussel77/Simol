#include "simol/statphys/dynamics/Langevin.hpp"

namespace simol
{
  ///Constructor for the Langevin Dynamics with thermostats everywhere
  Langevin::Langevin(Input const& input):
    LangevinBase(input)
  {
    //galerkin_ = createLangevinGalerkin(input);
    //galerkin_ = createGalerkin(input);
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
  void Langevin::computeGeneratorOnBasis(shared_ptr<CVBasis> cvBasis, System const& syst) const
  {
    cvBasis->generatorOnBasisValues_ = DVec::Zero(cvBasis->totalNbOfElts());
    //cout << "p : " << endl << cvBasis->pVariable(syst) << endl;
    //cout << "f : " << endl << cvBasis->forces(syst) << endl;
    for (int iOfFunction = 0; iOfFunction < cvBasis->totalNbOfElts(); iOfFunction++)
      //for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      cvBasis->generatorOnBasisValues_(iOfFunction) += dot(cvBasis->pVariable(syst) , cvBasis->gradientQ(syst, iOfFunction))
                              + dot(cvBasis->forces(syst) , cvBasis->gradientP(syst, iOfFunction))
                              + gamma() * (- dot(cvBasis->pVariable(syst) , cvBasis->gradientP(syst, iOfFunction))
                                        + cvBasis->laplacianP(syst, iOfFunction) / beta() );

      /*cout << "--iOfFunction : " << iOfFunction << endl;                        
      cout << "--gradQ : " << endl << cvBasis->gradientQ(syst, iOfFunction) << endl;
      cout << "--gradP : " << endl << cvBasis->gradientP(syst, iOfFunction) << endl;
      cout << "--laplaP : " << endl << cvBasis->laplacianP(syst, iOfFunction) << endl;*/
    }
    /*ofstream geneOnBasis("output/Langevin/Isolated/DoubleWell/HermiteHermite/geneOnBasis", std::ofstream::app);
    geneOnBasis << syst(0).position(0) << " " << syst(0).momentum(0) << " " << cvBasis->generatorOnBasisValues_(0) << " " << cvBasis->basisValues_[0] << endl;
  
    cout << "||||||||||||||OK---------------------" <<endl << endl;*/
  }
  
  
  
  ///Constructor for the Constrained Langevin Dynamics with thermostats everywhere
  ConstrainedLangevin::ConstrainedLangevin(Input const& input):
    Langevin(input),
    drift_(input.drift())
  {
    cout << "drift_ = " << drift_ << endl;
  }
  
  /*void ConstrainedLangevin::halfUpdateOrsteinUhlenbeck(System& syst) const
  {
    // Becareful here we assume that all the particles share the same mass !
    double alpha = exp(- gamma() / syst(0).mass() * timeStep()/2);
    double localLagrangeMultiplier = 0;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      localLagrangeMultiplier += 
      syst(iOfParticle).momentum() = alpha * particle.momentum() + sqrt((1 - pow(alpha, 2)) / localBeta * particle.mass()) * rng_->gaussian();
  }*/
  


}
