#include "tools.hpp"

double modulo(double variable, double mini, double maxi)
{
	double var2 = variable - mini;
	double maxi2 = maxi - mini;
	int ratio = floor(var2/maxi2);
	//cout << variable << " -> " << var2 - ratio * maxi2 + mini << endl;
	return var2 - ratio * maxi2 + mini;
}