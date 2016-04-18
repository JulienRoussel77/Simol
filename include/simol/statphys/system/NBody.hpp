#ifndef SIMOL_NBODY_HPP
#define SIMOL_NBODY_HPP
#include "System.hpp"

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
      list<int> const& members() const;
      list<int>& members();

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
      void computeAllForces();
      int nbOfParticlesPerDimension() const;
      double latticeParameter() const;
      void interaction(Particle& particle1, Particle& particle2) const;

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
