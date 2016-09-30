#include "simol/statphys/system/NBody.hpp"

namespace simol
{
  
  //--------------- particle pair iterators ------------
  
  
  
  /*const Cell& ParticlePairIterator::cell1() const 
  {return syst_->cell(iOfCell1());}
  
  Cell& ParticlePairIterator::cell1() 
  {return syst_->cell(iOfCell1());}  
  
  const int& ParticlePairIterator::iOfCell2() const 
  {return syst_->cell(iOfCell1()).indexNeighbors(iOfNeighbor2());}
  
  int& ParticlePairIterator::iOfCell2() 
  {return syst_->cell(iOfCell1()).indexNeighbors(iOfNeighbor2());}
  
  const Cell& ParticlePairIterator::cell2() const 
  {return syst_->cell(iOfCell2());}
  
  Cell& ParticlePairIterator::cell2() 
  {return syst_->cell(iOfCell2());}*/
  
  /*Particle const& ParticlePairIterator::particle2() const 
  {return syst_->getMember(iOfCell2(), iOfParticle2());}
  
  Particle& ParticlePairIterator::particle2() 
  {return syst_->getMember(iOfCell2(), iOfParticle2());}*/
  
  ParticlePairIterator NBody::pairBegin()
  {
    ParticlePairIterator iter = ParticlePairIterator();
    iter.it1_ = cell(0).members().begin();
    iter.it2_ = std::next(iter.it1_, 1);
    iter.endIt2_ = cell(0).members().end();
    return iter;
  }
  

  bool NBody::pairFinished(ParticlePairIterator const& it) const
  {
    return (it.iOfCell1() == nbOfCells());
  }
  
  void NBody::incrementePairIterator(ParticlePairIterator& it)
  {
    //if (cell(cell(it.iOfCell1()).indexNeighbors(it.iOfNeighbor2())).members().end() != it.endIt2()) cout << "endIt2 is not well defined !" << endl;
    it.it2_++;
    //if the neighbor cell is fully visited
    //if (it.it2_ == it.cell2().members().end())
    while (it.it2_ == it.endIt2())
    {
      //cout << "Box 2 finished !" << endl;
      //we jump to the next neighbor cell
      it.iOfNeighbor2()++;
      //Cell* cell2 = &(cell(cell(it.iOfCell1()).indexNeighbors(it.iOfNeighbor2())));
        
      if (it.iOfNeighbor2() != nbOfNeighbors())
        it.it2_ = getCell2(it).members().begin();
      //if all the neighbor cells has been visited
      else
      {
        //cout << "Loop over neighbor boxes finished !" << endl;
        //we jump to the next particle
        it.it1_++;
        it.iOfNeighbor2() = 0;
        //it.endIt2() = cell2->members().end();
        //if the cell is fully visited
        while (it.it1_ == cell(it.iOfCell1()).members().end())
        {
          //cout << "Box 1 finished !" << endl;
          //we jump to the next one
          it.iOfCell1()++;
          if (it.iOfCell1() != nbOfCells())
            it.it1_ = cell(it.iOfCell1()).members().begin();
          else
            return;
          //cout << "iOfCell1 / nbOfCells : " << it.iOfCell1() << " / " << nbOfCells() << endl;
        }
        //in both case this allows to avoid counting twice a pair or a particle with it self
        it.it2_ = std::next(it.it1_, 1);        
      }
      it.endIt2() = getCell2(it).members().end();
    }
    //if(!*it.it1_) cout << "it1_ is nullptr !" << endl;
    //if(!*it.it2_) cout << "it2_ is nullptr !" << endl;
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
    domainSize_(input.latticeParameter()*input.nbOfParticlesPerDimension()),
    restart_(input.restart()),
    doCells_(input.doCellMethod()),
    Rcut_(input.cutOffRatio()*input.potentialSigma())
  {
    //-- initialise the configuration by initializing the particles --
    assert(configuration_.size() > 1);
    for (int i = 0; i < input.nbOfParticles(); i++)
      configuration_[i] = new Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));

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

  bool NBody::restart() const
  {
    return restart_;
  }
  
  string NBody::restartFileName() const
  {
    return restartFileName_;
  }
  
  /*Particle& NBody::getMember(const int& iOfCell, const int& iOfMember)
  {return cell(iOfCell).member(iOfMember);}
  
  Particle const& NBody::getMember(const int& iOfCell, const int& iOfMember) const
  {return cell(iOfCell).member(iOfMember);}*/

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
      return returnIndexCell(i1, i2);
    else if (dimension_ == 3)
    {
      int i3 = floor(periodicPosition(pos(2)) / cellSize_);
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
    if (i2 != -1) indexCell += intModulo(i2, nbOfCellsPerDimension_) * nbOfCellsPerDimension_;
    if (i3 != -1) indexCell += intModulo(i3, nbOfCellsPerDimension_) * pow(nbOfCellsPerDimension_, 2);
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
    Vector<int> currentVecIndex(dimension_); 
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
  
  Vector<double> NBody::representant(Vector<double> vecDistance) const
  {
    Vector<double> periodicDistance;
    for (int d = 0; d < (int)dimension_; d++)
      periodicDistance(d) -= rint(vecDistance(d) / domainSize_) * domainSize_;
    return periodicDistance;
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

  ///
  /// compute the forces (with or without cell method)
  void NBody::computeAllForces()
  {
    //for (auto && particle : configuration_)
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).resetForce(potential());
    if (doCells_)
    {
      //-- reinitialize cells before looping on the pair interactions --
      reinitializeCells();
      //-- compute the interactions --
      for (ParticlePairIterator it = pairBegin(); !pairFinished(it); incrementePairIterator(it))
      {
        //cout << it.iOfCell1_ << " " << it.iOfNeighbor2_ << endl;
        interaction(it.particle1(), it.particle2());
      }
      
      /*int neighborIndex = 1;
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
      }*/
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





}
