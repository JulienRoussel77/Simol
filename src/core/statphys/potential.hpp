#ifndef SIMOL_POTENTIAL_HPP
#define SIMOL_POTENTIAL_HPP

#include "Vector.hpp"	
#include "input.hpp"

namespace simol
{
  class Potential;
  
  Potential* createPotential(Input const& input);

  class Potential
  {
    friend Potential* createPotential(Input const& input);
    
    public:
      Potential();
      virtual ~Potential(){};
      Potential(Input const& input);
      virtual double operator()(dvec const & position) const {std::cout << "This potential is pairwise !" << std::endl; exit(1);};
      virtual double operator()(double const & position) const {std::cout << "This potential is not pairwise !" << std::endl; exit(1);}
      virtual dvec derivative(dvec const & position) const = 0;
      virtual dvec force(dvec const & position) const = 0;
  };
  
  class Sinusoidal : public Potential{
    public:
      Sinusoidal(Input const& input);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;
      dvec force(dvec const & position) const;

    private:
      double amplitude_;
      double pulsation_;
  };
  
  class DoubleWell : public Potential{
    public:
      DoubleWell(Input const& input);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;
      dvec force(dvec const & position) const;

    private:
      double height_;
      double interWell_;
  };
  
    class Harmonic : public Potential{
    public:
      Harmonic(Input const& input);
      double operator()(double const & position) const;
      dvec derivative(dvec const & position) const;
      dvec force(dvec const & position) const;

    private:
      double stiffness_;
  };

}

//#include "potential.ipp"

#endif
