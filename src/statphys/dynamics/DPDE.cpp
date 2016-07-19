#include "simol/statphys/dynamics/DPDE.hpp"

#include "simol/statphys/dynamics/LangevinBase.hpp"

using std::cout;
using std::endl;

namespace simol
{
  DPDE::DPDE(Input const&  input):
    LangevinBase(input),
    heatCapacity_(input.heatCapacity()),
    cutOff_(input.cutOffRatio()*input.potentialSigma())
  {}

  void DPDE::printName() const
  {
    std::cout << "DynamicsType = DPDE" << std::endl;
  }

  const double& DPDE::heatCapacity() const
  {
    return heatCapacity_;
  }

  double& DPDE::heatCapacity()
  {
    return heatCapacity_;
  }

  double& DPDE::cutOff()
  {
    return cutOff_;
  }

  const double& DPDE::rejectionRate() const
  {
    return rejectionRate_;
  }
  
  double& DPDE::rejectionRate()
  {
    return rejectionRate_;
  }

  const double& DPDE::rejectionCount() const
  {
    return rejectionCount_;
  }
  
  double& DPDE::rejectionCount()
  {
    return rejectionCount_;
  }

  const double& DPDE::totalCountForRejection() const
  {
    return totalCountForRejection_;
  }
  
  double& DPDE::totalCountForRejection()
  {
    return totalCountForRejection_;
  }

  double DPDE::gamma_DPDE(double intEnergy)
  {
    return gamma_ * heatCapacity() * temperature() / intEnergy;
  }

  double DPDE::internalTemperature(double intEnergy) const
  {
    return intEnergy/heatCapacity();
  }

  double DPDE::sigma() const
  {
    return sqrt(2 * gamma_ * temperature());
  }

  void DPDE::incrementRejection()
  {
    rejectionCount_ += 1;
  }

  void DPDE::incrementTotalCountForRejection()
  {
    totalCountForRejection_ += 1;
  }

  //----------------- THERMALIZATION: LANGEVIN-LIKE FLUCT/DISS ---------------------

  //-- second part of Verlet + integrate the O-U analytically + update internal energies --
  void DPDE::secondPartThermalization(Particle& particle)
  {
    // second part of Verlet
    particle.momentum() += timeStep_ * particle.force() / 2;
    // Ornstein-Uhlenbeck process on the momenta
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1 - pow(alpha, 2)) / beta_ * particle.mass()) * rng_->gaussian();
    // update of the internal energies 
    particle.internalEnergy() += -(1-temperature()/internalTemperature(particle.internalEnergy()))*heatCapacity()*timeStep_ 
      + sqrt(2*temperature()*heatCapacity()*timeStep_)*rng_->scalarGaussian();
  }

  //-------------------------------- FOR ISOLATED SYSTEMS --------------------------
  
  void DPDE::energyReinjection(Particle& particle)
  {
    double old_kin_energy = particle.kineticEnergy();
    double local_gamma_DPDE = gamma_DPDE(particle.internalEnergy());
    double alpha = exp(- local_gamma_DPDE / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() 
      + sqrt((1 - pow(alpha, 2)) * particle.mass() * temperature() * gamma() / local_gamma_DPDE) * rng_->gaussian();
    double new_kin_energy = particle.kineticEnergy();
    particle.internalEnergy() += old_kin_energy - new_kin_energy;
  }

  Vector<double> DPDE::effectiveDrift(Particle& particle)
  {
    double e0 = particle.totalEnergyDPDE(); 
    Vector<double> drift = -pow(sigma(),2)* heatCapacity()/(2*particle.mass())*particle.momentum()/(e0-particle.kineticEnergy());
    return drift;
  }

  void DPDE::metropolizedEnergyReinjection(Particle& particle)
  {
    //-- keep previous configuration --
    Vector<double> old_momentum = particle.momentum();
    double old_kin_energy = particle.kineticEnergy();
    double E0 = old_kin_energy + particle.internalEnergy(); 
    //-- propose a new move --
    Vector<double> G = rng_->gaussian();
    // EM scheme: not stable...
    //Vector<double> b = effectiveDrift(particle);
    //particle.momentum() += b*timeStep_ + sigma()*sqrt(timeStep_)*G;
    // SSA like scheme: should be better
    double gamma_n = gamma_DPDE(E0-old_kin_energy);
    double alpha = exp(-gamma_n / particle.mass() * timeStep_);
    double sigma_n = sqrt((1 - pow(alpha, 2)) * particle.mass() * temperature() * gamma() / gamma_n);
    particle.momentum() = alpha * old_momentum + sigma_n * G;
    //-- compute the Metropolis rate --
    double new_kin_energy = particle.kineticEnergy();
    double rate = 0; 
    if (new_kin_energy < E0)
      {
	rate = heatCapacity()*(log(E0-new_kin_energy)-log(E0-old_kin_energy)) + 0.5*pow(G.norm(), 2); 
    	double gamma_reverse = gamma_DPDE(E0-new_kin_energy);
    	double alpha_reverse = exp(-gamma_reverse / particle.mass() * timeStep_);
    	double sigma2_reverse = (1 - pow(alpha_reverse,2)) * particle.mass() * temperature() * gamma() / gamma_reverse; 
    	Vector<double> GG = old_momentum - alpha_reverse*particle.momentum();
	rate -= 0.5*pow(GG.norm(), 2)/sigma2_reverse;
    	rate = exp(rate)*sigma_n/sqrt(sigma2_reverse);
      }
    //-- acceptance/rejection procedure --
    incrementTotalCountForRejection();
    double U = rng_->scalarUniform();
    if (U > rate)
      {
    	//-- reject the move --
    	incrementRejection();
	particle.momentum() = old_momentum;
      }
    else 
      {
	//-- update internal energies --
	particle.internalEnergy() = E0 - new_kin_energy;
      }
  }

  //-------------------------------- FOR ISOLATED SYSTEMS --------------------------
  
  double DPDE::chi(double dist)
  {
    return 1-dist/cutOff_;
  }

  double DPDE::pairwiseFluctuationDissipation(double v12, double dist, double internalEnergy1, double internalEnergy2, double reduced_mass)
  {
    // if the distance is beyond the cut off, no update
    if (dist > cutOff_)
      return v12;
    // otherwise, integrate the Ornstein-Uhlenbeck process of the projected velocity
    else 
      {
	double chi_n = chi(dist);
	double gamma_n = (gamma_DPDE(internalEnergy1) + gamma_DPDE(internalEnergy2)) / 2;
	double alpha_n = exp(- gamma_n * pow(chi_n,2)/reduced_mass * timeStep_);
	double G = rng_->scalarGaussian();
	double sigma_n = sqrt((1 - pow(alpha_n, 2)) * temperature() * gamma() / (gamma_n * reduced_mass ) );
	double new_v12 = alpha_n * v12 +  sigma_n * G;
	// update the rate for accept/reject correction
	rejectionRate() = 0.5*pow(G, 2) + log(sigma_n);
	// return the new relative velocity
	return new_v12;
      }
  }

  void DPDE::acceptRejectRate(double v12_current, double v12_init, double internalEnergy1, double internalEnergy2, double reduced_mass, double dist)
  {
    // Note that the rate is initialized in the function 'pairwiseFluctuationDissipation'
    //cout << min(internalEnergy1,internalEnergy2) << " : " << internalEnergy1 << " and " << internalEnergy2 << endl;
    double Delta_v2 = (pow(v12_current,2) - pow(v12_init,2))*reduced_mass/4;
    if (Delta_v2 < min(internalEnergy1,internalEnergy2))
      {
	// rejectionRate() = 1; // always accept -- for debugging
	// compute the new energies
	double E1 = internalEnergy1 - Delta_v2;
	double E2 = internalEnergy2 - Delta_v2;
	rejectionRate() += heatCapacity()*(log(E1) + log(E2) - log(internalEnergy1) - log(internalEnergy2));
	double chi_n = chi(dist);
	double gamma_reverse = ( gamma_DPDE(E1) + gamma_DPDE(E2) ) / 2;
	double alpha_reverse = exp(- gamma_reverse * pow(chi_n,2)/reduced_mass * timeStep_);
	double sigma_reverse = sqrt((1 - pow(alpha_reverse, 2)) * temperature() * gamma() / (gamma_reverse*reduced_mass) );
	double GG = (alpha_reverse*v12_current - v12_init)/sigma_reverse;
	rejectionRate() -= 0.5*pow(GG, 2) + log(sigma_reverse);
	rejectionRate() = exp(rejectionRate());
      }
    //-- otherwise systematically reject --
    else 
      rejectionRate() = 0;
  }

}
