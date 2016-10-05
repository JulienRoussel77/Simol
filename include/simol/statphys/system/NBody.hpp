#ifndef SIMOL_NBODY_HPP
#define SIMOL_NBODY_HPP
#include "simol/statphys/system/System.hpp"

namespace simol
{

  //----- when the cell method is used -----
  class Cell
  {
    
  public:
    
    Cell();
    int size() const;
    void clear();
    void push_back(Particle& particle);
    vector<int> const& indexNeighbors() const;
    vector<int>& indexNeighbors();
    int const& indexNeighbors(const int& iOfNeighbor) const;
    int& indexNeighbors(const int& iOfNeighbor);
    vector<Particle*> const& members() const;
    vector<Particle*>& members();
    Particle const& member(const int& iOfMember) const;
    Particle& member(const int& iOfMember);
    int nbOfMembers() const;
    
  protected:
    
    vector<Particle*> members_;
    vector<int> indexNeighbors_;
    
  };
  
  //------------------ class for multiparticle systems in 2D/3D ------------------------
  class NBody : public System
  {

  public:

    NBody(Input const& input);
    void printName() const;
    
    //-- pair iterator --
    ParticlePairIterator pairBegin();
    bool pairFinished(ParticlePairIterator const& it) const;
    void incrementePairIterator(ParticlePairIterator& it);
    Cell const& getCell2(ParticlePairIterator const& it) const;
    Cell& getCell2(ParticlePairIterator const& it);
    
    // -- accessors --
    int const& nbOfNeighbors() const {return nbOfNeighbors_;}
    int const& nbOfCells() const {return nbOfCells_;}
    Cell const& cell(int const& iOfCell) const {return cells_[iOfCell];}
    Cell & cell(int const& iOfCell) {return cells_[iOfCell];}
  
    virtual Vector<double> representant(Vector<double> vecDistance) const;
    void samplePositions(DynamicsParameters const& dynaPara);
    void computeAllForces();
    int nbOfParticlesPerDimension() const;
    double latticeParameter() const;
    void interaction(Particle& particle1, Particle& particle2) const;
  
    //--- for cell method ---
    void reinitializeCells();
    void findNeighboringCells();
    int findIndex(Vector<double> const& pos) const;
    int returnIndexCell(int i1, int i2 = -1, int i3 = -1) const;
    double periodicPosition(double x) const;

  //protected:
  public:
    int nbOfParticlesPerDimension_;
    double latticeParameter_;
    double domainSize_;

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
