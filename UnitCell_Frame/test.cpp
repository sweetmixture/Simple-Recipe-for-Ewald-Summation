#include <iostream>
#include "cell.hpp"

#define TARGET_INPUT "geo.txt"

int main()
{
	Cell c(TARGET_INPUT);
	c.CalcCoulombEnergy();
	c.CalcCoulombDerivative();
	c.ShowBasicCellInfo();
	return 0;
}
