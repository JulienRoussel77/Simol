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

bool hasSmallerNorm(const cplx& a, const cplx& b)
{
  return (norm(a) < norm(b));
}