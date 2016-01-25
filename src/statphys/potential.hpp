#ifndef SIMOL_POTENTIAL_HPP
#define SIMOL_POTENTIAL_HPP

#include "tools.hpp"
#include "core/linalg/Vector.hpp"
#include "input.hpp"


namespace simol
{
  class Potential;

  Potential* createPotential(Input const& input, int const& indexOfReplica=1);

  class Potential
  {
    friend Potential* createPotential(Input const& input, int const& indexOfReplica);

    public:
      Potential();
      virtual ~Potential(){};
      Potential(Input const& input, int const& indexOfReplica=1);
      virtual double operator()(dvec const & /*position*/) const {std::cout << "This potential is pairwise !" << std::endl; exit(1);};
      virtual dvec derivative(dvec const & position) const = 0;
      virtual double laplacian(dvec const & /*position*/) const {std::cout << "Laplacian not implemented !" << std::endl; exit(1);};
      virtual dvec force(dvec const & position) const;
  };

  class Sinusoidal : public Potential{
    public:
      Sinusoidal(Input const& input, int const& indexOfReplica=1);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;
      double laplacian(dvec const & position) const;

    private:
      double amplitude_;
      double pulsation_;
  };

    class SumSinusoidal : public Potential{
    public:
      SumSinusoidal(Input const& input, int const& indexOfReplica=1);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;
      virtual double laplacian(dvec const & position) const;


    private:
      double amplitude_;
      double pulsation_;
  };

  class FracSinusoidal : public Potential{
    public:
      FracSinusoidal(Input const& input, int const& indexOfReplica=1);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;
      virtual double laplacian(dvec const & position) const;


    private:
      double amplitude_;
      double pulsation_;
  };

  class DoubleWell : public Potential{
    public:
      DoubleWell(Input const& input, int const& indexOfReplica=1);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;

    private:
      double height_;
      double interWell_;
  };

    class HarmonicWell : public Potential{
    public:
      HarmonicWell(Input const& input, int const& indexOfReplica=1);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;

    private:
      double stiffness_;
  };




  class Harmonic : public Potential{
    public:
      Harmonic(Input const& input, int const& indexOfReplica=1);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;

    private:
      double stiffness_;
  };

  class Rotor : public Potential{
    public:
      Rotor(Input const& input, int const& indexOfReplica=1);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;
      double laplacian(dvec const & position) const;

    private:
      double stiffness_;
  };


    class Quadratic : public Potential{
    public:
      Quadratic(Input const& input, int const& indexOfReplica=1);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;

    private:
      double stiffness_, alpha_, beta_;
  };
}

//#include "potential.ipp"

#endif
