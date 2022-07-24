#ifndef __INTERACTION_Manager
#define __INTERACTION_Manager

#include <iostream>

#include "atom.hpp"

class Cell;

class Manager // Interaction Manager
{

public:
	void info( const Cell& cell )
	{
		std::cout << "Manager\n";
		std::cout << "CellVol : " << cell.volume << std::endl;
	}
	
	double mono_mono_real( const Atom& atom_i, const Atom& atom_j )
	{
		double res = 0.;
		return res;
	}







};

#endif
