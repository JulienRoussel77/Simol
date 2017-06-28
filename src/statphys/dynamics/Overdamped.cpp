#include "simol/statphys/dynamics/Overdamped.hpp"

namespace simol
{
  //#### Overdamped ####

  ///
  ///Constructor for the Overdamped Langevin Dynamics
  Overdamped::Overdamped(Input const& input):
    Dynamics(input)
  {
    //galerkin_ = createOverdampedGalerkin(input);
    //galerkin_ = createGalerkin(input);
  }
 

  // /!\ a changer dans le cas N != 1
  ///Proceeds to an Euler-Maruyama scheme (in comment the second order is available)
  void Overdamped::updatePosition(Particle& particle)
  {
    //particle.position() += timeStep_ * particle.force() + sqrt(2 * timeStep_ / beta_) * rng_->gaussian();
    DVec newGaussian = rng_->gaussian();
    particle.position() += timeStep_ * particle.force() + sqrt(timeStep_ / (2*beta_)) * (newGaussian + particle.oldGaussian());
    particle.oldGaussian() = newGaussian;
  }
  
  void Overdamped::getThermo(Output& output) const
  {    
    double kineticEnergy = output.dimension() * output.nbOfParticles() / (2 * beta());
    
    output.temperature() = 1 / beta();
    output.totalEnergy() = kineticEnergy + output.potentialEnergy();
  }
  
  ///
  ///Computes the pressure from the kineticEnergy and the totalVirial, these fields must be updated
  void Overdamped::getPressure(Output& output) const
  {
    double kineticEnergy = output.dimension() * output.nbOfParticles() / (2 * beta());
    output.pressure() = (2 * kineticEnergy + output.totalVirial()) / (output.dimension() * output.nbOfParticles() * pow(output.latticeParameter(), output.dimension()));
  }

  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  void Overdamped::computeGeneratorOnBasis(shared_ptr<CVBasis> cvBasis, System const& syst) const
  {
    /*cout << "computeGeneratorOnBasis : totalNbOfElts = " << cvBasis->totalNbOfElts() << endl;
    //cout << "p : " << endl << cvBasis->pVariable(syst) << endl;
    cout << "r : " << endl << cvBasis->qVariable(syst) << endl;
    cout << "f : " << endl << cvBasis->forces(syst) << endl;*/
    cvBasis->generatorOnBasisValues_ = DVec::Zero(cvBasis->totalNbOfElts());
    //DVec result = DVec::Zero(nbOfFunctions());
    DMat forceMat = cvBasis->forces(syst);
    for (int iOfFunction = 0; iOfFunction < cvBasis->totalNbOfElts(); iOfFunction++)
    {
      //for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      cvBasis->generatorOnBasisValues_(iOfFunction) += cvBasis->laplacianQ(syst, iOfFunction) / beta()
                               + dot(forceMat, cvBasis->gradientQ(syst, iOfFunction));
      //cvBasis->generatorOnBasisValues_(iOfFunction) += (cvBasis->forces(syst))(0,0);                        
      //cvBasis->generatorOnBasisValues_(iOfFunction) += (cvBasis->gradientQ(syst, iOfFunction))(0,0);
      /*cout << "--iOfFunction : " << iOfFunction << endl;
      cout << "--laplacianQ : " << cvBasis->laplacianQ(syst, iOfFunction) << endl;
      cout << "--gradientQ : " << cvBasis->gradientQ(syst, iOfFunction) << endl;*/
    }
    //ofstream test("testest", std::ofstream::app);
    //test << cvBasis->qVariable(syst) << " " << forceMat << " " << cvBasis->basisValues_[0] << " " << cvBasis->gradientQ(syst, 0) << " " << cvBasis->laplacianQ(syst, 0) << endl;
                   
    //ofstream tempOut("generatorOnBasis.txt", std::ofstream::app);
    //tempOut << syst(0).position(0) << " " << syst(0).force(0) << " " << dot(*cvBasis->cvCoeffs_, cvBasis->generatorOnBasisValues_) << endl;
  }

}
