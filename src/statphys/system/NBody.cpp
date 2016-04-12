#include "NBody.hpp"

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

  // ------------------- NBody ----------------------------
  
  NBody::NBody(Input const& input):
    System(input),
    nbOfParticlesPerDimension_(input.nbOfParticlesPerDimension()),
    latticeParameter_(input.latticeParameter()),
    domainSize_(input.latticeParameter()*input.nbOfParticlesPerDimension()),
    Rcut_(input.cutOffRatio()*input.potentialSigma())
  {
    //-- initialise the configuration by initializing the particles --
    assert(configuration_.size() > 1);
    for (int i = 0; i<input.nbOfParticles(); i++) 
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
    
    //-- initialize the cells --
    nbOfCellsPerDimension_ = floor(domainSize_/Rcut_ );
    assert(nbOfCellsPerDimension_ >= 3);
    nbOfCells_ = pow(nbOfCellsPerDimension_,dimension_);
    cells_ = vector<Cell>(nbOfCells_, Cell());
    //cout << " size: " << domainSize_ << ", cut off: " <<  input.cutOffRatio()*input.potentialSigma() << ", nb cells per dim = " << nbOfCellsPerDimension_ << ", total = " << nbOfCells_ << endl;
    
    if (dimension_ == 2)
      nbOfNeighbors_ = 4;
    if (dimension_ == 3)
      nbOfNeighbors_ = 13;
    
    findNeighboringCells();
    reinitializeCells();
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

  ///
  /// return position x in the interval [0,domainSize_)
  double NBody::periodicPosition(double x) const
  {
    return x - floor(x/domainSize_)*domainSize_;
  }
  
  int NBody::findIndex(Vector<double> const& pos) const
  {
    // PBM 2 : appel a pos(0)
    int i1 = floor(periodicPosition(pos(0))/Rcut_);
    int i2 = floor(periodicPosition(pos(1))/Rcut_);
    if (dimension_ == 2)
      return returnIndexCell2D(i1,i2);
    else if (dimension_ == 3)
      {
    	int i3 = floor(periodicPosition(pos(2))/Rcut_);
    	return returnIndexCell3D(i1,i2,i3);
      }
    else
      throw std::invalid_argument("NBody::findIndex -- dimension incorrect");
    return 0;
  }

  int NBody::returnIndexCell2D(int i1, int i2) const
  {
    //cout << "i1 = " << i1 << ", modulo = " << intModulo(i1,nbOfCellsPerDimension_) 
    //	 << "; i2 = " << i2 << ", modulo = " << intModulo(i2,nbOfCellsPerDimension_) << endl;
    return  intModulo(i1,nbOfCellsPerDimension_) + intModulo(i2,nbOfCellsPerDimension_) * nbOfCellsPerDimension_;
  }

  int NBody::returnIndexCell3D(int i1, int i2, int i3) const
  {
    return  intModulo(i1,nbOfCellsPerDimension_) + nbOfCellsPerDimension_ * intModulo(i2,nbOfCellsPerDimension_) 
      + pow(nbOfCellsPerDimension_,2) * intModulo(i3,nbOfCellsPerDimension_); 
  }
  
  void NBody::findNeighboringCells()
  {
    //-- initialize the vector of indices of neighboring cells --
    for (int j = 0; j < nbOfCells_; j++)
      cells_[j].indexNeighbors() = vector<int>(nbOfNeighbors_,0);
    //-- find the neighboring cells --
    int currentCellIndex;
    if (dimension_ == 2)
      {
	for (int i1 = 0; i1 < nbOfCellsPerDimension_; i1++)
	  {
	    for (int i2 = 0; i2 < nbOfCellsPerDimension_; i2++)
	      {
		currentCellIndex = returnIndexCell2D(i1,i2);
		//cout << i1 << " " << i2 << " : " << currentCellIndex << endl;
		cells_[currentCellIndex].indexNeighbors()[0] = returnIndexCell2D(i1+1,i2);
		cells_[currentCellIndex].indexNeighbors()[1] = returnIndexCell2D(i1+1,i2+1);
		cells_[currentCellIndex].indexNeighbors()[2] = returnIndexCell2D(i1,  i2+1);
		cells_[currentCellIndex].indexNeighbors()[3] = returnIndexCell2D(i1-1,i2+1);
		// cout << "  ";
		// for (int k = 0; k < 4; k++)
		//   cout << cells_[currentCellIndex].indexNeighbors()[k] << " ";
		// cout << endl;
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
		     currentCellIndex = returnIndexCell3D(i1,i2,i3);
		     cout << i1 << " " << i2 << " " << i3 << " : " << currentCellIndex << endl;
		     //-- 4 vertices out of 8 --
		     cells_[currentCellIndex].indexNeighbors()[0]  = returnIndexCell3D(i1+1,i2+1,i3+1);
		     cells_[currentCellIndex].indexNeighbors()[1]  = returnIndexCell3D(i1+1,i2+1,i3-1);
		     cells_[currentCellIndex].indexNeighbors()[2]  = returnIndexCell3D(i1+1,i2-1,i3+1);
		     cells_[currentCellIndex].indexNeighbors()[3]  = returnIndexCell3D(i1+1,i2-1,i3-1);
		     //-- 3 faces out of 6 --
		     cells_[currentCellIndex].indexNeighbors()[4]  = returnIndexCell3D(i1+1,i2  ,i3);
		     cells_[currentCellIndex].indexNeighbors()[5]  = returnIndexCell3D(i1  ,i2+1,i3);
		     cells_[currentCellIndex].indexNeighbors()[6]  = returnIndexCell3D(i1  ,i2  ,i3+1);
		     //-- 6 edges out of 12 --
		     cells_[currentCellIndex].indexNeighbors()[7]  = returnIndexCell3D(i1  ,i2+1,i3+1);
		     cells_[currentCellIndex].indexNeighbors()[8]  = returnIndexCell3D(i1+1,i2  ,i3+1);
		     cells_[currentCellIndex].indexNeighbors()[9]  = returnIndexCell3D(i1+1,i2+1,i3);
		     cells_[currentCellIndex].indexNeighbors()[10] = returnIndexCell3D(i1  ,i2+1,i3-1);
		     cells_[currentCellIndex].indexNeighbors()[11] = returnIndexCell3D(i1+1,i2  ,i3-1);
		     cells_[currentCellIndex].indexNeighbors()[12] = returnIndexCell3D(i1+1,i2-1,i3);
		     cout << "  ";
		     for (int k = 0; k < 13; k++)
		       cout << cells_[currentCellIndex].indexNeighbors()[k] << " ";
		     cout << endl;
		   }
	      }
	  }
      }
  }

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

  void NBody::computeAllForces()
  {
    for (auto&& particle : configuration_)
      particle.resetForce(potential());
    //-- reinitialize cells before looping on the pair interactions --
    //reinitializeCells()
    for (int i = 0; i < nbOfParticles(); i++)
      for (int j = i+1; j < nbOfParticles(); j++)
      interaction(configuration_[i], configuration_[j]);
  }
    
  void NBody::interaction(Particle& particle1, Particle& particle2) const
  {
    Vector<double> r12 = particle1.position() - particle2.position();
    // take closest periodic image
    double distance = 0.;
    for (int d = 0; d < (int)dimension_; d++)
      {
        r12(d) -= rint(r12(d)/domainSize_)*domainSize_;
        distance += r12(d)*r12(d);
      }
    distance = sqrt(distance);
    // compute energy
    double energy12 = potential(distance);
    // cout << distance << "  " << energy12 << endl;
    particle1.potentialEnergy() += energy12/2;
    particle2.potentialEnergy() += energy12/2;
    // compute forces 
    double force12 = potentialForce(distance)(0);
    r12 /= distance;
    particle1.force() += force12 * r12;
    particle2.force() -= force12 * r12;
    // compute pressure, based on the Virial formula for the potential part: P_pot = -\sum_{i < j} r_ij v'(r_ij) / d|Vol|; will divide by d|Vol| at the end
    particle1.virial() += 0.5*force12*distance;
    particle2.virial() += 0.5*force12*distance;
  }
  
}
