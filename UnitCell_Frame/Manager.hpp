#ifndef __INTERACTION_Manager
#define __INTERACTION_Manager

#include <iostream>

//#include "atom.hpp"

class Cell;
class Atom;

class Manager // Interaction Manager
{

public:
	void info( const Cell& cell )
	{
		std::cout << "Manager\n";
		std::cout << "CellVol : " << cell.volume << std::endl;
	}
	
	double coulomb_mono_mono_real( const Cell& cell, const Atom& atom_i, const Atom& atom_j )
	{
		double q_prod = atom_i.charge*atom_j.charge;
			


// this->mono_real_energy += 0.5*(this->AtomList[i]->charge*this->AtomList[j]->charge)/r_norm * erfc(r_norm/this->sigma) * this->TO_EV;

		std::cout << atom_i.type << std::endl;
		std::cout << atom_j.type << std::endl;

		double res = 0.;
		return res;
	}

};

#endif
