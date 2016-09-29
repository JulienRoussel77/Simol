#ifndef SIMOL_NBODY_HPP
#define SIMOL_NBODY_HPP
#include "simol/statphys/system/System.hpp"
#include "simol/statphys/dynamics/DPDE.hpp"

namespace simol
{

  //----- when the cell method is used -----
  class Cell
  {
    
  public:
    
    Cell();
    int size() const;
    void clear();
    void push_back(int ind);
    vector<int> const& indexNeighbors() const;
    vector<int>& indexNeighbors();
    int const& indexNeighbors(const int& iOfNeighbor) const;
    int& indexNeighbors(const int& iOfNeighbor);
    list<int> const& members() const;
    list<int>& members();
    int nbOfMembers() const;
    
  protected:
    
    list<int> members_;
    vector<int> indexNeighbors_;
    
  };
  
  //------------------ class for multiparticle systems in 2D/3D ------------------------
  class NBody : public System
  {

  public:

    NBody(Input const& input);
    void printName() const;
    
    //-- pair iterator --
    bool pairFinished(ParticlePairIterator const& it) const;
    void incrementePairIterator(ParticlePairIterator& it);
    
    // -- accessors --
    int const& nbOfNeighbors() const {return nbOfNeighbors_;}
    int const& nbOfCells() const {return nbOfCells_;}
    Cell const& cell(int const& iOfCell) const {return cells_[iOfCell];}
    Cell & cell(int const& iOfCell) {return cells_[iOfCell];}
  
    void computeAllForces();
    void fluctuationDissipationDPDE(DPDE& dyna);
    int nbOfParticlesPerDimension() const;
    double latticeParameter() const;
    void interaction(Particle& particle1, Particle& particle2) const;
    void elementaryFluctuationDissipationDPDE(DPDE& dyna, Particle& particle1, Particle& particle2) const;

    //--- for initial conditions (restart) --
    bool restart() const;
    string restartFileName() const;
  
    //--- for cell method ---
    void reinitializeCells();
    void findNeighboringCells();
    int findIndex(Vector<double> const& pos) const;
    int returnIndexCell2D(int i1, int i2) const;
    int returnIndexCell3D(int i1, int i2, int i3) const;
    double periodicPosition(double x) const;

  protected:

    int nbOfParticlesPerDimension_;
    double latticeParameter_;
    double domainSize_;

    //--- if restarting from a given file ---
    bool restart_;
    string restartFileName_;

    //--- for cell method ---
    bool doCells_;
    double Rcut_;
    double cellSize_;
    int nbOfCellsPerDimension_;
    int nbOfCells_;
    int nbOfNeighbors_;
    vector<Cell> cells_;

  };

}


#endif
