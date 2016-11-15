#include "simol/statphys/dynamics/DPDE.hpp"

#include "simol/statphys/dynamics/LangevinBase.hpp"

using std::cout;
using std::endl;

namespace simol
{
  DPDE::DPDE(Input const&  input):
    LangevinBase(input),
    heatCapacity_(input.heatCapacity()),
    einsteinTemperature_(input.einsteinTemperature()),
    kappa_(input.kappa()),
    cutOff_(input.cutOffRatio()*input.potentialSigma()),
    rejectionCountFD_(0),
    totalCountForRejectionFD_(0),
    negativeEnergiesCountFD_(0),
    rejectionCountThermal_(0),
    totalCountForRejectionThermal_(0),
    negativeEnergiesCountThermal_(0),
    rejectionRate_(0),
    doMetropolis_(input.doMetropolis())
  {}

  const double& DPDE::heatCapacity() const
  {
    return heatCapacity_;
  }

  double& DPDE::heatCapacity()
  {
    return heatCapacity_;
  }

  const double& DPDE::einsteinTemperature() const
  {
    return einsteinTemperature_;
  }

  double& DPDE::einsteinTemperature()
  {
    return einsteinTemperature_;
  }

  const double& DPDE::kappa() const
  {
    return kappa_;
  }

  double& DPDE::kappa()
  {
    return kappa_;
  }

  double& DPDE::cutOff()
  {
    return cutOff_;
  }

  //-- for 1D systems --
  double DPDE::gamma_DPDE(double intEnergy)
  {
    return gamma_ * heatCapacity() * temperature() / intEnergy;
  }

  double DPDE::sigma() const
  {
    return sqrt(2 * gamma_ * temperature());
  }

  //-- for Metropolis procedures --
  const double& DPDE::rejectionRate() const {return rejectionRate_;}
  double& DPDE::rejectionRate() {return rejectionRate_;}
  // FD part
  const double& DPDE::rejectionCountFD() const {return rejectionCountFD_;}
  double& DPDE::rejectionCountFD() {return rejectionCountFD_;}
  const double& DPDE::negativeEnergiesCountFD() const {return negativeEnergiesCountFD_;}
  double& DPDE::negativeEnergiesCountFD() {return negativeEnergiesCountFD_;}
  const double& DPDE::totalCountForRejectionFD() const {return totalCountForRejectionFD_;}
  double& DPDE::totalCountForRejectionFD() {return totalCountForRejectionFD_;}
  void DPDE::incrementRejectionFD() {rejectionCountFD_ += 1;}
  void DPDE::incrementTotalCountForRejectionFD() {totalCountForRejectionFD_ += 1;}
  // Thermal part 
  const double& DPDE::rejectionCountThermal() const {return rejectionCountThermal_;}
  double& DPDE::rejectionCountThermal() {return rejectionCountThermal_;}
  const double& DPDE::negativeEnergiesCountThermal() const {return negativeEnergiesCountThermal_;}
  double& DPDE::negativeEnergiesCountThermal() {return negativeEnergiesCountThermal_;}
  const double& DPDE::totalCountForRejectionThermal() const {return totalCountForRejectionThermal_;}
  double& DPDE::totalCountForRejectionThermal() {return totalCountForRejectionThermal_;}
  void DPDE::incrementRejectionThermal() {rejectionCountThermal_ += 1;}
  void DPDE::incrementTotalCountForRejectionThermal() {totalCountForRejectionThermal_ += 1;}

  //-------------------- micro EOS for NBody systems ----------------------------
  
  double DPDE::entropy(double intEnergy) const
  {
    if (einsteinTemperature() > 0)
      return ((intEnergy+heatCapacity()*einsteinTemperature())*log(intEnergy+heatCapacity()*einsteinTemperature()) - intEnergy*log(intEnergy) )/einsteinTemperature(); 
    else
      return heatCapacity()*log(intEnergy);
  }

  double DPDE::entropy_derivative(double intEnergy) const
  {
    if (einsteinTemperature() > 0)
      return -log( intEnergy/(intEnergy+heatCapacity()*einsteinTemperature()) )/einsteinTemperature();
    else
      return heatCapacity()/intEnergy;
  }

  double DPDE::internalTemperature(double intEnergy) const
  {
    return 1./entropy_derivative(intEnergy);
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
    //---- update of the internal energies; requires metropolization; also use timestep "adapted" to cv. rate of the dynamics ---- 
    double effectiveTimeStep = timeStep_*heatCapacity();
    double G = rng_->scalarGaussian();
    double old_energy = particle.internalEnergy();
    particle.internalEnergy() += -(1-entropy_derivative(particle.internalEnergy())*temperature())*effectiveTimeStep
      + sqrt(2*temperature()*effectiveTimeStep)*G;
    double rate = 0; 
    double GG = 0;
    if (particle.internalEnergy() > 0)
      {
	rate += entropy(particle.internalEnergy()) - entropy(old_energy) - (particle.internalEnergy() - old_energy)/temperature();
	GG = (old_energy - particle.internalEnergy() + (1-entropy_derivative(particle.internalEnergy())*temperature())*effectiveTimeStep)/sqrt(2*temperature()*effectiveTimeStep);
	rate += 0.5*( pow(G,2) - pow(GG,2) );
	rate = exp(rate);
      }
    double U = rng_->scalarUniform();
    if (U > rate) //-- reject the move --
      particle.internalEnergy() = old_energy;
  }

  //-------------------------------- FLUCT/DISS FOR ISOLATED SYSTEMS --------------------------
  
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

  DVec DPDE::effectiveDrift(Particle& particle)
  {
    double e0 = particle.totalEnergyDPDE(); 
    DVec drift = -pow(sigma(),2)* heatCapacity()/(2*particle.mass())*particle.momentum()/(e0-particle.kineticEnergy());
    return drift;
  }

  void DPDE::metropolizedEnergyReinjection(Particle& particle)
  {
    //-- keep previous configuration --
    DVec old_momentum = particle.momentum();
    double old_kin_energy = particle.kineticEnergy();
    double E0 = old_kin_energy + particle.internalEnergy(); 
    //-- propose a new move --
    DVec G = rng_->gaussian();
    // EM scheme: not stable...
    //DVec b = effectiveDrift(particle);
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
    	DVec GG = old_momentum - alpha_reverse*particle.momentum();
	rate -= 0.5*pow(GG.norm(), 2)/sigma2_reverse;
    	rate = exp(rate)*sigma_n/sqrt(sigma2_reverse);
      }
    //-- acceptance/rejection procedure --
    incrementTotalCountForRejectionFD();
    double U = rng_->scalarUniform();
    if (U > rate)
      {
    	//-- reject the move --
    	incrementRejectionFD();
	particle.momentum() = old_momentum;
      }
    else 
      {
	//-- update internal energies --
	particle.internalEnergy() = E0 - new_kin_energy;
      }
  }

  //-------------------------------- FLUCT/DISS FOR NBODY SYSTEMS --------------------------
  
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
	double gamma_n = 0.5 * gamma() * temperature() * (entropy_derivative(internalEnergy1) + entropy_derivative(internalEnergy2));
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
	rejectionRate() += entropy(E1) + entropy(E2) - entropy(internalEnergy1) - entropy(internalEnergy2); 
	double chi_n = chi(dist);
	double gamma_reverse = 0.5 * gamma() * temperature() * (entropy_derivative(E1) + entropy_derivative(E2));
	double alpha_reverse = exp(- gamma_reverse * pow(chi_n,2)/reduced_mass * timeStep_);
	double sigma_reverse = sqrt((1 - pow(alpha_reverse, 2)) * temperature() * gamma() / (gamma_reverse*reduced_mass) );
	double GG = (alpha_reverse*v12_current - v12_init)/sigma_reverse;
	rejectionRate() -= 0.5*pow(GG, 2) + log(sigma_reverse);
	rejectionRate() = exp(rejectionRate());
      }
    //-- otherwise systematically reject --
    else 
      {
	rejectionRate() = 0;
	negativeEnergiesCountFD() += 1;
      }
  }

  ///
  /// compute the fluctuation/dissipation in DPDE using SSA
  void DPDE::fluctuationDissipation(System& syst)
  {
    if (syst.doCells())
      {
	//-- reinitialize cells before looping on the pair interactions --
	syst.reinitializeCells();
	//-- compute the interactions --
	for (ParticlePairIterator it = syst.pairBegin(); !syst.pairFinished(it); syst.incrementePairIterator(it))
	  elementaryFluctuationDissipation(syst, it.particle1(), it.particle2());
      }
    else
      {
	//-- no cell method: std double loop --
	for (int i = 0; i < syst.nbOfParticles(); i++)
	  for (int j = i + 1; j < syst.nbOfParticles(); j++)
	    elementaryFluctuationDissipation(syst, syst(i), syst(j));
      }
  }
  
  ///
  /// elementary interaction between two particles
  //----------- !!!! UPDATE FORMULAS FOR PARTICLES WITH DIFFERENT MASSES !!!! ------------
  void DPDE::elementaryFluctuationDissipation(System const& syst, Particle& particle1, Particle& particle2)
  {
    // keep previous configuration
    DVec old_momentum_1 = particle1.momentum();
    DVec old_momentum_2 = particle2.momentum();
    // compute the unit vector e12 of line of centers and the distance (as above) 
    DVec r12 = syst.periodicImage(particle1.position() - particle2.position());
    double distance = r12.norm();
    // if the distance is not below the cutoff, nothing has to be done
    if (distance < cutOff_)
      {
	DVec e12 = r12/distance;
	// compute the variation of the relative velocity
	double mu12 = 1./( 1./particle1.mass() + 1./particle2.mass() ); // reduced mass
	DVec vect12 = particle1.momentum()/particle1.mass() - particle2.momentum()/particle2.mass();  
	double v12_0 = dot(vect12,e12);
	double v12 = pairwiseFluctuationDissipation(v12_0,distance,particle1.internalEnergy(),particle2.internalEnergy(),mu12); 
	// update the momenta
	DVec totalMomentum = particle1.momentum() + particle2.momentum();
	DVec v12_perp = vect12 - dot(vect12,e12)*e12;
	particle1.momentum() += mu12*(v12-v12_0)*e12;
	particle2.momentum() -= mu12*(v12-v12_0)*e12;
	//mu12*( totalMomentum/particle2.mass() + v12_perp + v12*e12);
	//particle2.momentum() = mu12*( totalMomentum/particle1.mass() - v12_perp - v12*e12);
	// accept/reject step
	incrementTotalCountForRejectionFD();
	acceptRejectRate(v12,v12_0,particle1.internalEnergy(),particle2.internalEnergy(),mu12,distance);
	// check whether Metropolis required
	double U = -1; // in order to accept even if energies are negative...
	if (doMetropolis_)
	  U = rng_->scalarUniform();
	if (U > rejectionRate())
	  {
	    //-- reject the move --
	    incrementRejectionFD();
	    particle1.momentum() = old_momentum_1;
	    particle2.momentum() = old_momentum_2;
	  }
	else 
	  {
	    //-- update internal energies --
	    double internal_energy_variation = mu12*( pow(v12,2)-pow(v12_0,2) )/4;
	    particle1.internalEnergy() -= internal_energy_variation;
	    particle2.internalEnergy() -= internal_energy_variation;
	  }
      }
  }
  
  //--------------------------- THERMAL CONDUCTION ------------------------------

  ///
  /// compute the thermal conduction in DPDE using SSA
  void DPDE::thermalConduction(System& syst)
  {
    if (syst.doCells())
      {
	//-- reinitialize cells before looping on the pair interactions --
	syst.reinitializeCells();
	//-- compute the interactions --
	for (ParticlePairIterator it = syst.pairBegin(); !syst.pairFinished(it); syst.incrementePairIterator(it))
	  elementaryThermalConduction(syst, it.particle1(), it.particle2());
      }
    else
      {
	//-- no cell method: std double loop --
	for (int i = 0; i < syst.nbOfParticles(); i++)
	  for (int j = i + 1; j < syst.nbOfParticles(); j++)
	    elementaryThermalConduction(syst, syst(i), syst(j));
      }
  }

  ///
  /// elementary thermal conduction exchange between two particles
  void DPDE::elementaryThermalConduction(System const& syst, Particle& particle1, Particle& particle2)
  {
    double old_internalEnergy_1 = particle1.internalEnergy();
    double old_internalEnergy_2 = particle2.internalEnergy();
    // compute the distance 
    DVec r12 = syst.periodicImage(particle1.position() - particle2.position());
    double distance = r12.norm();
    if (distance < cutOff_)
      {
    	// compute predicted energy variation
    	double chi_n = chi(distance);
    	double deltaInternalEnergy = pairwiseThermalConduction(chi_n,old_internalEnergy_1,old_internalEnergy_2);
	// check whether the energy variation is admissible; otherwise set acceptance rate to 0
    	if ( (old_internalEnergy_1+deltaInternalEnergy>0) & (old_internalEnergy_2-deltaInternalEnergy>0)) 
    	  {
    	    particle1.internalEnergy() += deltaInternalEnergy;
    	    particle2.internalEnergy() -= deltaInternalEnergy;
    	    rejectionRate() += entropy(particle1.internalEnergy())+entropy(particle2.internalEnergy())-entropy(old_internalEnergy_1)-entropy(old_internalEnergy_2);
	    // compute the probability of the reverse move
    	    double GG = old_internalEnergy_1 - particle1.internalEnergy()
    	      - kappa()*pow(chi_n,2)*(entropy_derivative(particle1.internalEnergy())-entropy_derivative(particle2.internalEnergy()))*timeStep_ ;
    	    GG /= chi_n*sqrt(2*kappa()*timeStep_);
    	    rejectionRate() += -0.5*pow(GG,2);
	    rejectionRate() = exp(rejectionRate());
    	  }
    	else
    	  {
    	    rejectionRate() = 0;
    	    negativeEnergiesCountThermal() += 1;
    	  }
	// actual acceptance/rejection step
    	incrementTotalCountForRejectionThermal();
	double U = -1;  // in order to accept even if energies are negative...
	if (doMetropolis_)
	  U = rng_->scalarUniform();
	if (U > rejectionRate())
    	  {
	    //-- reject the move --
    	    incrementRejectionThermal();
    	    particle1.internalEnergy() = old_internalEnergy_1;
    	    particle2.internalEnergy() = old_internalEnergy_2;
    	  }
      }
  }

  double DPDE::pairwiseThermalConduction(double chi_n, double internalEnergy1, double internalEnergy2)
  {
    double G = rng_->scalarGaussian();
    double energy_increment = kappa()*pow(chi_n,2)*(entropy_derivative(internalEnergy1)-entropy_derivative(internalEnergy2))*timeStep_ 
      + sqrt(2*kappa()*timeStep_)*chi_n*G;
    // update the rate for accept/reject correction
    rejectionRate() = 0.5*pow(G, 2);
    // return the energy increment
    return energy_increment;
  }

  //------------------------ Specific output functions ------------------------------------

  void DPDE::getThermo(Output& output) const
  {
    output.temperature() = 2 * output.kineticEnergy() / (output.dimension() * output.nbOfParticles());
    output.totalEnergy() = output.kineticEnergy() + output.potentialEnergy() + output.internalEnergy();
  }
  
  void DPDE::specificComputeOutput(Output& output) const
  {
    //-- rejection rates --
    output.rejectionCountFD() = rejectionCountFD()/totalCountForRejectionFD();
    output.negativeEnergiesCountFD() = negativeEnergiesCountFD();
    output.rejectionCountThermal() = rejectionCountThermal()/totalCountForRejectionThermal();
    output.negativeEnergiesCountThermal() = negativeEnergiesCountThermal();
  }
  
}
