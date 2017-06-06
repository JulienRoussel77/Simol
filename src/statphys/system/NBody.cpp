#include "simol/statphys/system/NBody.hpp"

namespace simol
{
  
  //--------------- particle pair iterators ------------
  
  
  /*ParticlePairIterator NBody::pairBegin()
  {
    if (nbOfParticles() < 2) throw runtime_error("Cannot initialize a Particle Pair Iterator if the system has 0 or 1 particle!");
    ParticlePairIterator iter = ParticlePairIterator();
    int iOfCell = 0;
    iter.it1_ = cell(0).members().begin();
    // Look for the first non empty cell
    while (iter.it1_ == cell(iOfCell).members().end())
    {
      iOfCell++;
      iter.it1_ = cell(iOfCell).members().begin();
    }
    iOfCell1_ = iOfCell;
    iter.it2_ = std::next(iter.it1_, 1);
    int iOfNeighbour = 1;
    // If the first cell has a single particle, look for the next non empty cell
    while (iter.it2_ == cell(iOfCell).members().end())
    {
      iOfCell++;
      iter.it2_ = cell(iOfCell).members().begin();
    }
    iter.endIt2_ = cell(iOfCell).members().end();
    return iter;
  }*/
  
  ParticlePairIterator NBody::pairBegin()
  {
    if (nbOfParticles() < 2) throw runtime_error("Cannot initialize a Particle Pair Iterator if the system has 0 or 1 particle!");
    ParticlePairIterator iter = ParticlePairIterator();
    int iOfCell = 0;
    iter.it1_ = cell(0).members().begin();
    // Look for the first non empty cell
    while (iter.it1_ == cell(iOfCell).members().end())
    {
      iOfCell++;
      iter.it1_ = cell(iOfCell).members().begin();
    }
    iter.iOfCell1_ = iOfCell;
    iter.it2_ = iter.it1_;
    iter.iOfNeighbor2() = 0;
    iter.endIt2() = cell(iOfCell).members().end();
    
    incrementePairIterator(iter);
    return iter;
  }
  
  bool NBody::pairFinished(ParticlePairIterator const& it) const
  {
    return (it.iOfCell1() == nbOfCells());
  }
  
  void NBody::incrementePairIterator(ParticlePairIterator& it)
  {
    it.it2_++;
    //if the neighbor cell is fully visited
    while (it.it2_ == it.endIt2())
    {
      //we jump to the next neighbor cell
      it.iOfNeighbor2()++;
        
      if (it.iOfNeighbor2() != nbOfNeighbors())
        it.it2_ = getCell2(it).members().begin();
      //if all the neighboring cells of particle 1 have been visited, increment iterator of particle 1
      else
      {
        //we jump to the next particle 1
        it.it1_++;
        it.iOfNeighbor2() = 0;
        //if the cell is fully visited
        while (it.it1_ == cell(it.iOfCell1()).members().end())
        {
          //we jump to the next one
          it.iOfCell1()++;
          if (it.iOfCell1() != nbOfCells())
            it.it1_ = cell(it.iOfCell1()).members().begin();
          else
            return;
        }
        //in both case this allows to avoid counting twice a pair or a particle with it self
        it.it2_ = std::next(it.it1_, 1);        
      }
      it.endIt2() = getCell2(it).members().end();
    }
  }
  
  Cell const& NBody::getCell2(ParticlePairIterator const& it) const
  {
    return cell(cell(it.iOfCell1()).indexNeighbors(it.iOfNeighbor2()));
  }
  
  Cell& NBody::getCell2(ParticlePairIterator const& it)
  {
    return cell(cell(it.iOfCell1()).indexNeighbors(it.iOfNeighbor2()));
  }
  
  
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

  void Cell::push_back(Particle& particle)
  {
    members_.push_back(&particle);
  }

  vector<int> const& Cell::indexNeighbors() const
  {
    return indexNeighbors_;
  }

  vector<int>& Cell::indexNeighbors()
  {
    return indexNeighbors_;
  }
  
  int const& Cell::indexNeighbors(const int& iOfNeighbor) const
  {return indexNeighbors_[iOfNeighbor];}
  
  int& Cell::indexNeighbors(const int& iOfNeighbor)
  {return indexNeighbors_[iOfNeighbor];}
  
  vector<Particle*> const& Cell::members() const
  {
    return members_;
  }

  vector<Particle*>& Cell::members()
  {
    return members_;
  }
  
  Particle const& Cell::member(int const& iOfMember) const
  {
    return *(members_[iOfMember]);
  }

  Particle& Cell::member(int const& iOfMember)
  {
    return *(members_[iOfMember]);
  }
  
  int Cell::nbOfMembers() const
  {
    return members_.size();
  }

  // ------------------- NBody ----------------------------

  NBody::NBody(Input const& input):
    System(input),
    nbOfParticlesPerDimension_(input.nbOfParticlesPerDimension()),
    latticeParameter_(input.latticeParameter()),
    doCells_(input.doCellMethod()),
    Rcut_(input.cutOffRatio()*input.potentialSigma())
  {
    domainSize_ = input.latticeParameter()*input.nbOfParticlesPerDimension();
    
    //-- initialise the configuration by initializing the particles, intial__ reads the restart file if it exists --
    assert(configuration_.size() > 1);
    for (int iOfParticle = 0; iOfParticle < input.nbOfParticles(); iOfParticle++)
      configuration_[iOfParticle] = new Particle(input.mass(), input.initialPosition(iOfParticle), input.initialMomentum(iOfParticle), input.initialInternalEnergy(iOfParticle));

    //-- initialize the cells --
    if (doCells_)
    {
      // creation of the cells (only if more than 3 per direction)
      nbOfCellsPerDimension_ = floor(domainSize_ / Rcut_ );
      if (nbOfCellsPerDimension_ < 3)
        throw runtime_error("There must be at least 3 cells per direction ! Disable the cells and rerun.");
      cellSize_ = domainSize_ / nbOfCellsPerDimension_;
      nbOfCells_ = pow(nbOfCellsPerDimension_, dimension_);
      cells_ = vector<Cell>(nbOfCells_, Cell());
      // number of neighboring cells: depends on the dimension
      if (dimension_ == 2)
        nbOfNeighbors_ = 4+1;
      if (dimension_ == 3)
        nbOfNeighbors_ = 13+1;
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

  ///
  /// find the index of the cell in which a particle is
  int NBody::findIndex(DVec const& pos) const
  {
    DVec perPos = periodicImage(pos);
    int i1 = floor(perPos(0) / cellSize_);
    int i2 = floor(perPos(1) / cellSize_);
    if (dimension_ == 2)
      return returnIndexCell(i1, i2);
    else if (dimension_ == 3)
    {
      int i3 = floor(perPos(2) / cellSize_);
      return returnIndexCell(i1, i2, i3);
    }
    else
      throw std::invalid_argument("NBody::findIndex -- dimension incorrect");
    return 0;
  }
  
  ///
  /// determine the index of a cell from the vector of integers characterizing the cell
  int NBody::returnIndexCell(int i1, int i2, int i3) const
  {
              int indexCell  = intModulo(i1, nbOfCellsPerDimension_);
    if (i2 != -1) indexCell = indexCell * nbOfCellsPerDimension_ + intModulo(i2, nbOfCellsPerDimension_);
    if (i3 != -1) indexCell = indexCell * nbOfCellsPerDimension_ + intModulo(i3, nbOfCellsPerDimension_);
    return indexCell;
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
    if (dimension_ == 2)
    {
      for (int i1 = 0; i1 < nbOfCellsPerDimension_; i1++)
      {
        for (int i2 = 0; i2 < nbOfCellsPerDimension_; i2++)
        {
          int currentCellIndex = returnIndexCell(i1, i2);
          cells_[currentCellIndex].indexNeighbors(0) = returnIndexCell(i1    , i2);
          cells_[currentCellIndex].indexNeighbors(1) = returnIndexCell(i1 + 1, i2);
          cells_[currentCellIndex].indexNeighbors(2) = returnIndexCell(i1 + 1, i2 + 1);
          cells_[currentCellIndex].indexNeighbors(3) = returnIndexCell(i1,  i2 + 1);
          cells_[currentCellIndex].indexNeighbors(4) = returnIndexCell(i1 - 1, i2 + 1);
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
            int currentCellIndex = returnIndexCell(i1, i2, i3);
            //-- current cell itself --
            cells_[currentCellIndex].indexNeighbors(0)  = returnIndexCell(i1      , i2    , i3    );
            //-- 4 vertices out of 8 --
            cells_[currentCellIndex].indexNeighbors(1)  = returnIndexCell(i1 + 1, i2 + 1, i3 + 1);
            cells_[currentCellIndex].indexNeighbors(2)  = returnIndexCell(i1 + 1, i2 + 1, i3 - 1);
            cells_[currentCellIndex].indexNeighbors(3)  = returnIndexCell(i1 + 1, i2 - 1, i3 + 1);
            cells_[currentCellIndex].indexNeighbors(4)  = returnIndexCell(i1 + 1, i2 - 1, i3 - 1);
            //-- 3 faces out of 6 --
            cells_[currentCellIndex].indexNeighbors(5)  = returnIndexCell(i1 + 1, i2  , i3);
            cells_[currentCellIndex].indexNeighbors(6)  = returnIndexCell(i1  , i2 + 1, i3);
            cells_[currentCellIndex].indexNeighbors(7)  = returnIndexCell(i1  , i2  , i3 + 1);
            //-- 6 edges out of 12 --
            cells_[currentCellIndex].indexNeighbors(8)  = returnIndexCell(i1  , i2 + 1, i3 + 1);
            cells_[currentCellIndex].indexNeighbors(9)  = returnIndexCell(i1 + 1, i2  , i3 + 1);
            cells_[currentCellIndex].indexNeighbors(10)  = returnIndexCell(i1 + 1, i2 + 1, i3);
            cells_[currentCellIndex].indexNeighbors(11) = returnIndexCell(i1  , i2 + 1, i3 - 1);
            cells_[currentCellIndex].indexNeighbors(12) = returnIndexCell(i1 + 1, i2  , i3 - 1);
            cells_[currentCellIndex].indexNeighbors(13) = returnIndexCell(i1 + 1, i2 - 1, i3);
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
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    {
      int cellIndex = findIndex(getParticle(iOfParticle).position());
      cell(cellIndex).push_back(getParticle(iOfParticle));
    }
  }
  
  
  
  void NBody::samplePositions(DynamicsParameters const& dynaPara)
  {
    cout << " - Sampling the positions..." << endl;
    int Dim = dimension();
    int NbPartDim = nbOfParticlesPerDimension();
    double latticeSize = latticeParameter();

    for (int i = 0; i < nbOfParticles(); i++)
    {
      getParticle(i).momentum() = drawMomentum(dynaPara.beta(), getParticle(i).mass());
      getParticle(i).countdown() = i;
    }
    //-- initialization on a cubic lattice --
    if (Dim == 2)
    {
      for (int i = 0; i < NbPartDim; i++)
        for (int j = 0; j < NbPartDim; j++)
        {
          getParticle(i * NbPartDim + j).position(0) = (i+.5) * latticeSize;
          getParticle(i * NbPartDim + j).position(1) = (j+.5) * latticeSize;
        }
    }
    else if (Dim == 3)
    {
      int NbPartDim2 = NbPartDim * NbPartDim;
      for (int i = 0; i < NbPartDim; i++)
        for (int j = 0; j < NbPartDim; j++)
          for (int k = 0; k < NbPartDim; k++)
          {
            getParticle(i * NbPartDim2 + j * NbPartDim + k).position(0) = (i+.5) * latticeSize;
            getParticle(i * NbPartDim2 + j * NbPartDim + k).position(1) = (j+.5) * latticeSize;
            getParticle(i * NbPartDim2 + j * NbPartDim + k).position(2) = (k+.5) * latticeSize;
          }
    }
    else
    {
      throw std::invalid_argument("sampleSystem: Bad dimension, should be 2 or 3");
    }
  }
  

  ///
  /// compute the forces (with or without cell method)
  void NBody::computeAllForces()
  {
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).resetForce(externalPotential());
    
    if (doCells_)
    {
      //-- reinitialize cells before looping on the pair interactions --
      reinitializeCells();
      //-- compute the interactions --
      for (ParticlePairIterator it = pairBegin(); !pairFinished(it); incrementePairIterator(it))
        interaction(it.particle1(), it.particle2());
    }
    else
    {
      //-- no cell method: std double loop --
      for (int i = 0; i < nbOfParticles(); i++)
        for (int j = i + 1; j < nbOfParticles(); j++)
          interaction(getParticle(i), getParticle(j));
    }
  }

  ///
  /// elementary interaction between two particles
  void NBody::interaction(Particle& particle1, Particle& particle2) const
  {
    // take closest periodic image
    /*cout << "####################" << endl;
    cout << particle1.position().adjoint() << " <-> " << particle2.position().adjoint() << endl;*/
    DVec r12 = periodicDistance(particle1.position() - particle2.position());
    //cout << r12.adjoint() << endl;
    double distance = r12.norm();
    /*cout << distance << endl;
    cout << particle1.force().adjoint() << endl;
    cout << "####################" << endl;*/

    // compute energy
    double energy12 = pairPotential()(distance);
    particle1.potentialEnergy() += energy12 / 2;
    particle2.potentialEnergy() += energy12 / 2;
    // compute forces
    double force12 = pairPotential().scalarPotentialForce(distance);
    r12 /= distance;
    particle1.force() += force12 * r12;
    particle2.force() -= force12 * r12;
    // compute pressure, based on the Virial formula for the potential part: P_pot = -\sum_{i < j} r_ij v'(r_ij) / d|Vol|; will divide by d|Vol| at the end
    particle1.virial() += 0.5 * force12 * distance;
    particle2.virial() += 0.5 * force12 * distance;
    //cout << distance << " -> " << force12 << " " << energy12 << endl;

    //ofstream pairs("pairs.txt", std::ofstream::app);
    //pairs << distance << endl;
    
  }





}
