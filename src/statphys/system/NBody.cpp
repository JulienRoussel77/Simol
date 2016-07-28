#include "simol/statphys/system/NBody.hpp"

namespace simol
{
  //-------------------- cells ----------------------

  Cell::Cell()
  {
  }

  int Cell::size() const
  {
    return members_.size();
  }

  void Cell::clear()
  {
    members_.clear();
  }

  void Cell::push_back(int ind)
  {
    members_.push_back(ind);
  }

  vector<int> const& Cell::indexNeighbors() const
  {
    return indexNeighbors_;
  }

  vector<int>& Cell::indexNeighbors()
  {
    return indexNeighbors_;
  }

  list<int> const& Cell::members() const
  {
    return members_;
  }

  list<int>& Cell::members()
  {
    return members_;
  }

  // ------------------- NBody ----------------------------

  NBody::NBody(Input const& input):
    System(input),
    nbOfParticlesPerDimension_(input.nbOfParticlesPerDimension()),
    latticeParameter_(input.latticeParameter()),
    domainSize_(input.latticeParameter()*input.nbOfParticlesPerDimension()),
    restart_(input.restart()),
    doCells_(input.doCellMethod()),
    Rcut_(input.cutOffRatio()*input.potentialSigma())
  {
    //-- initialise the configuration by initializing the particles --
    assert(configuration_.size() > 1);
    for (int i = 0; i < input.nbOfParticles(); i++)
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));

    //-- check whether the initial configuration should be re-read from a file --
    if (restart_)
      {
	restartFileName_ = input.restartFileName();
      }

    //-- initialize the cells --
    if (doCells_)
    {
      // creation of the cells (only if more than 3 per direction)
      nbOfCellsPerDimension_ = floor(domainSize_ / Rcut_ );
      assert(nbOfCellsPerDimension_ >= 3);
      cellSize_ = domainSize_ / nbOfCellsPerDimension_;
      nbOfCells_ = pow(nbOfCellsPerDimension_, dimension_);
      cells_ = vector<Cell>(nbOfCells_, Cell());
      // number of neighboring cells: depends on the dimension
      if (dimension_ == 2)
        nbOfNeighbors_ = 4;
      if (dimension_ == 3)
        nbOfNeighbors_ = 13;
      // create the array of neighboring cell indices
      findNeighboringCells();
      // create lists of particles in each cell
      reinitializeCells();
    }
  }

  void NBody::printName() const
  {
    cout << "SystemType = NBody" << endl;
  }

  int NBody::nbOfParticlesPerDimension() const
  {
    return nbOfParticlesPerDimension_;
  }

  double NBody::latticeParameter() const
  {
    return latticeParameter_;
  }

  bool NBody::restart() const
  {
    return restart_;
  }
  
  string NBody::restartFileName() const
  {
    return restartFileName_;
  }

  ///
  /// return position x in the interval [0,domainSize_)
  double NBody::periodicPosition(double x) const
  {
    return x - floor(x / domainSize_) * domainSize_;
  }

  ///
  /// find the index of the cell in which a particle is
  int NBody::findIndex(Vector<double> const& pos) const
  {
    int i1 = floor(periodicPosition(pos(0)) / cellSize_);
    int i2 = floor(periodicPosition(pos(1)) / cellSize_);
    if (dimension_ == 2)
      return returnIndexCell2D(i1, i2);
    else if (dimension_ == 3)
    {
      int i3 = floor(periodicPosition(pos(2)) / cellSize_);
      return returnIndexCell3D(i1, i2, i3);
    }
    else
      throw std::invalid_argument("NBody::findIndex -- dimension incorrect");
    return 0;
  }

  ///
  /// determine the index of a cell in a 2D setting from the vector of integers characterizing the cell
  int NBody::returnIndexCell2D(int i1, int i2) const
  {
    return  intModulo(i1, nbOfCellsPerDimension_) + intModulo(i2, nbOfCellsPerDimension_) * nbOfCellsPerDimension_;
  }

  ///
  /// determine the index of a cell in a 3D setting from the vector of integers characterizing the cell
  int NBody::returnIndexCell3D(int i1, int i2, int i3) const
  {
    return  intModulo(i1, nbOfCellsPerDimension_) + nbOfCellsPerDimension_ * intModulo(i2, nbOfCellsPerDimension_)
            + pow(nbOfCellsPerDimension_, 2) * intModulo(i3, nbOfCellsPerDimension_);
  }

  ///
  /// construct the list of neighboring cells for all cells
  /// only half of the neighbors are considered (in order to avoid double counting of interactions)
  void NBody::findNeighboringCells()
  {
    //-- initialize the vector of indices of neighboring cells --
    for (int j = 0; j < nbOfCells_; j++)
      cells_[j].indexNeighbors() = vector<int>(nbOfNeighbors_, 0);
    //-- find the neighboring cells --
    int currentCellIndex;
    if (dimension_ == 2)
    {
      for (int i1 = 0; i1 < nbOfCellsPerDimension_; i1++)
      {
        for (int i2 = 0; i2 < nbOfCellsPerDimension_; i2++)
        {
          currentCellIndex = returnIndexCell2D(i1, i2);
          cells_[currentCellIndex].indexNeighbors()[0] = returnIndexCell2D(i1 + 1, i2);
          cells_[currentCellIndex].indexNeighbors()[1] = returnIndexCell2D(i1 + 1, i2 + 1);
          cells_[currentCellIndex].indexNeighbors()[2] = returnIndexCell2D(i1,  i2 + 1);
          cells_[currentCellIndex].indexNeighbors()[3] = returnIndexCell2D(i1 - 1, i2 + 1);
        }
      }
    }
    if (dimension_ == 3)
    {
      for (int i1 = 0; i1 < nbOfCellsPerDimension_; i1++)
      {
        for (int i2 = 0; i2 < nbOfCellsPerDimension_; i2++)
        {
          for (int i3 = 0; i3 < nbOfCellsPerDimension_; i3++)
          {
            currentCellIndex = returnIndexCell3D(i1, i2, i3);
            //-- 4 vertices out of 8 --
            cells_[currentCellIndex].indexNeighbors()[0]  = returnIndexCell3D(i1 + 1, i2 + 1, i3 + 1);
            cells_[currentCellIndex].indexNeighbors()[1]  = returnIndexCell3D(i1 + 1, i2 + 1, i3 - 1);
            cells_[currentCellIndex].indexNeighbors()[2]  = returnIndexCell3D(i1 + 1, i2 - 1, i3 + 1);
            cells_[currentCellIndex].indexNeighbors()[3]  = returnIndexCell3D(i1 + 1, i2 - 1, i3 - 1);
            //-- 3 faces out of 6 --
            cells_[currentCellIndex].indexNeighbors()[4]  = returnIndexCell3D(i1 + 1, i2  , i3);
            cells_[currentCellIndex].indexNeighbors()[5]  = returnIndexCell3D(i1  , i2 + 1, i3);
            cells_[currentCellIndex].indexNeighbors()[6]  = returnIndexCell3D(i1  , i2  , i3 + 1);
            //-- 6 edges out of 12 --
            cells_[currentCellIndex].indexNeighbors()[7]  = returnIndexCell3D(i1  , i2 + 1, i3 + 1);
            cells_[currentCellIndex].indexNeighbors()[8]  = returnIndexCell3D(i1 + 1, i2  , i3 + 1);
            cells_[currentCellIndex].indexNeighbors()[9]  = returnIndexCell3D(i1 + 1, i2 + 1, i3);
            cells_[currentCellIndex].indexNeighbors()[10] = returnIndexCell3D(i1  , i2 + 1, i3 - 1);
            cells_[currentCellIndex].indexNeighbors()[11] = returnIndexCell3D(i1 + 1, i2  , i3 - 1);
            cells_[currentCellIndex].indexNeighbors()[12] = returnIndexCell3D(i1 + 1, i2 - 1, i3);
          }
        }
      }
    }
  }

  ///
  /// recreate the list of particles for each cell
  void NBody::reinitializeCells()
  {
    //-- clear the lists --
    for (int j = 0; j < nbOfCells_; j++)
      cells_[j].clear();
    //-- reaffect the particles --
    int index = 0;
    for (int i = 0; i < nbOfParticles(); i++)
    {
      index = findIndex(getParticle(i).position());
      cells_[index].push_back(i);
    }
  }

  ///
  /// compute the forces (with or without cell method)
  void NBody::computeAllForces()
  {
  for (auto && particle : configuration_)
      particle.resetForce(potential());
    if (doCells_)
    {
      //-- reinitialize cells before looping on the pair interactions --
      reinitializeCells();
      //-- compute the interactions --
      int neighborIndex = 1;
      for (int k = 0; k < nbOfCells_; k++)
      {
        //-- interaction within cells: avoid double counting by setting i1 \leq i2 + 1 --
        for (list<int>::iterator it1 = cells_[k].members().begin(); it1 != cells_[k].members().end(); it1++)
          for (list<int>::iterator it2 = std::next(it1, 1); it2 != cells_[k].members().end(); it2++)
            interaction(configuration_[*it1], configuration_[*it2]);
        //-- interactions between neighboring cells: full double loops --
        for (int l = 0; l < nbOfNeighbors_; l++)
        {
          // index of neighboring cell
          neighborIndex = cells_[k].indexNeighbors()[l];
          // complete double loop between the elements of cells_[k] and its neighbor cells_[neighborIndex]
          for (list<int>::iterator it1 = cells_[k].members().begin(); it1 != cells_[k].members().end(); it1++)
            for (list<int>::iterator it2 = cells_[neighborIndex].members().begin(); it2 != cells_[neighborIndex].members().end(); it2++)
              interaction(configuration_[*it1], configuration_[*it2]);
        }
      }
    }
    else
    {
      //-- no cell method: std double loop --
      for (int i = 0; i < nbOfParticles(); i++)
        for (int j = i + 1; j < nbOfParticles(); j++)
          interaction(configuration_[i], configuration_[j]);
    }
  }

  ///
  /// elementary interaction between two particles
  void NBody::interaction(Particle& particle1, Particle& particle2) const
  {
    Vector<double> r12 = particle1.position() - particle2.position();
    // take closest periodic image
    double distance = 0.;
    for (int d = 0; d < (int)dimension_; d++)
    {
      r12(d) -= rint(r12(d) / domainSize_) * domainSize_;
      distance += r12(d) * r12(d);
    }
    distance = sqrt(distance);
    // compute energy
    double energy12 = potential(distance);
    particle1.potentialEnergy() += energy12 / 2;
    particle2.potentialEnergy() += energy12 / 2;
    // compute forces
    double force12 = potentialForce(distance)(0);
    r12 /= distance;
    particle1.force() += force12 * r12;
    particle2.force() -= force12 * r12;
    // compute pressure, based on the Virial formula for the potential part: P_pot = -\sum_{i < j} r_ij v'(r_ij) / d|Vol|; will divide by d|Vol| at the end
    particle1.virial() += 0.5 * force12 * distance;
    particle2.virial() += 0.5 * force12 * distance;
  }

  ///
  /// compute the fluctuation/dissipation in DPDE using SSA
  void NBody::fluctuationDissipationDPDE(DPDE& dyna)
  {
    if (doCells_)
      {
	//-- reinitialize cells before looping on the pair interactions --
	reinitializeCells();
	//-- compute the interactions --
	int neighborIndex = 1;
	for (int k = 0; k < nbOfCells_; k++)
	  {
	    //-- interaction within cells: avoid double counting by setting i1 \leq i2 + 1 --
	    for (list<int>::iterator it1 = cells_[k].members().begin(); it1 != cells_[k].members().end(); it1++)
	      for (list<int>::iterator it2 = std::next(it1, 1); it2 != cells_[k].members().end(); it2++)
		elementaryFluctuationDissipationDPDE(dyna, configuration_[*it1], configuration_[*it2]);
	    //-- interactions between neighboring cells: full double loops --
	    for (int l = 0; l < nbOfNeighbors_; l++)
	      {
		// index of neighboring cell
		neighborIndex = cells_[k].indexNeighbors()[l];
		// complete double loop between the elements of cells_[k] and its neighbor cells_[neighborIndex]
		for (list<int>::iterator it1 = cells_[k].members().begin(); it1 != cells_[k].members().end(); it1++)
		  for (list<int>::iterator it2 = cells_[neighborIndex].members().begin(); it2 != cells_[neighborIndex].members().end(); it2++)
		    elementaryFluctuationDissipationDPDE(dyna, configuration_[*it1], configuration_[*it2]);
	      }
	  }
      }
    else
      {
	//-- no cell method: std double loop --
	for (int i = 0; i < nbOfParticles(); i++)
	  for (int j = i + 1; j < nbOfParticles(); j++)
	    elementaryFluctuationDissipationDPDE(dyna, configuration_[i], configuration_[j]);
      }
  }

  ///
  /// elementary interaction between two particles
  //----------- !!!! UPDATE FORMULAS FOR PARTICLES WITH DIFFERENT MASSES !!!! ------------
  void NBody::elementaryFluctuationDissipationDPDE(DPDE& dyna, Particle& particle1, Particle& particle2) const
  {
    // keep previous configuration
    Vector<double> old_momentum_1 = particle1.momentum();
    Vector<double> old_momentum_2 = particle2.momentum();
    // compute the unit vector e12 of line of centers and the distance (as above) 
    Vector<double> r12 = particle1.position() - particle2.position();
    double distance = 0.;
    for (int d = 0; d < (int)dimension_; d++)
    {
      r12(d) -= rint(r12(d) / domainSize_) * domainSize_;
      distance += r12(d) * r12(d);
    }
    distance = sqrt(distance);
    Vector<double> e12 = r12/distance;
    // compute the variation of the relative velocity
    //double old_kin_energy = particle1.kineticEnergy() + particle2.kineticEnergy();
    //dyna.printName(); 
    double mu12 = 1./( 1./particle1.mass() + 1./particle2.mass() ); // reduced mass
    Vector<double> vect12 = particle1.momentum()/particle1.mass() - particle2.momentum()/particle2.mass();  
    double v12_0 = dot(vect12,e12);
    double v12 = dyna.pairwiseFluctuationDissipation(v12_0,distance,particle1.internalEnergy(),particle2.internalEnergy(),mu12); 
    // update the momenta
    Vector<double> totalMomentum = particle1.momentum() + particle2.momentum();
    Vector<double> v12_perp = vect12 - dot(vect12,e12)*e12;
    //cout << vect12 << ", " << e12 << " : dot = " << dot(vect12,e12) << endl;
    particle1.momentum() += mu12*(v12-v12_0)*e12;
    particle2.momentum() -= mu12*(v12-v12_0)*e12;
    //mu12*( totalMomentum/particle2.mass() + v12_perp + v12*e12);
    //particle2.momentum() = mu12*( totalMomentum/particle1.mass() - v12_perp - v12*e12);
    // accept/reject step
    dyna.incrementTotalCountForRejection();
    dyna.acceptRejectRate(v12,v12_0,particle1.internalEnergy(),particle2.internalEnergy(),mu12,distance);
    double U = rng_->scalarUniform();
    if (U > dyna.rejectionRate())
      {
    	//-- reject the move --
    	dyna.incrementRejection();
	particle1.momentum() = old_momentum_1;
	particle2.momentum() = old_momentum_2;
      }
    else 
      {
	//-- update internal energies --
	//double new_kin_energy = particle1.kineticEnergy() + particle2.kineticEnergy();
	//double internal_energy_variation = 0.5*(new_kin_energy-old_kin_energy);
	double internal_energy_variation = mu12*( pow(v12,2)-pow(v12_0,2) )/4;
	//cout << mu12*( pow(v12,2)-pow(v12_0,2) )/4 - internal_energy_variation << endl;
	particle1.internalEnergy() -= internal_energy_variation;
	particle2.internalEnergy() -= internal_energy_variation;
      }
  }

}
