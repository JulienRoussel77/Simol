#include "simol/statphys/Tools.hpp"

double modulo(double variable, double mini, double maxi)
{
  double var2 = variable - mini;
  double maxi2 = maxi - mini;
  int ratio = floor(var2 / maxi2);
  //cout << variable << " -> " << var2 - ratio * maxi2 + mini << endl;
  return var2 - ratio * maxi2 + mini;
}

int intModulo(int variable, int maxi)
{
  return variable - floor(variable * 1. / maxi) * maxi;
}

void displayTime(double time)
{
  std::cout << endl << "End of simulation" << std::endl;
  int nbSeconds = time / CLOCKS_PER_SEC;
  int nbMinutes = nbSeconds / 60;
  int nbHours = nbMinutes / 60;
  int nbDays = nbHours / 24;
  cout << "Total computation time: " << nbDays << " j " <<  nbHours % 24  << " h " << nbMinutes % 60 << " m " << nbSeconds % 60 << " s " << endl;
}

int getNbOfLines(ifstream& file)
{
  std::string line;
  int nbOfLines = 0;
  for (nbOfLines = 0; std::getline(file, line); ++nbOfLines) {}
  return nbOfLines;
}

/*bool hasSmallerNorm(const cplx& a, const cplx& b)
{
  return (norm(a) < norm(b));
}*/

bool hasSmallerNorm(cplx a, cplx b)
{
  return (norm(a) < norm(b));
}

double dot(DVec const& u, DVec const& v)
{
  return u.dot(v);
}

DMat reshape(DVec const& u, int nbOfRows, int nbOfCols)
{
  DMat A(nbOfRows, nbOfCols);
  for (int j=0; j<nbOfCols; j++)
    for (int i=0; i<nbOfRows; i++)
      A(i,j) = u(j * nbOfRows + i);
  return A;
}

DVec rint(DVec const& u)
{
  DVec v(u.size());
  for (int d = 0; d < (int)u.rows(); d++)
    v(d) = rint(u(d));
  return v;
}

DVec polynomialProduct(DVec const& P, DVec const& Q)
{
  DVec PQ(DVec::Zero(P.rows() + Q.rows() - 1));
  for (int iOfCoeffP = 0; iOfCoeffP < P.rows(); iOfCoeffP++)
    for (int iOfCoeffQ = 0; iOfCoeffQ < Q.rows(); iOfCoeffQ++)
    {
      //cout << iOfCoeff << " = " << iOfCoeffP << " + " << iOfCoeff - iOfCoeffP << " -> " << P(iOfCoeffP) << " x " << Q(iOfCoeff - iOfCoeffP) << endl;
      PQ(iOfCoeffP + iOfCoeffQ) += P(iOfCoeffP) * Q(iOfCoeffQ);
    }
  return PQ;
}


