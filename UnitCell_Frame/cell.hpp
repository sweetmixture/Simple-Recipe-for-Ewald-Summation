#ifndef __CELL
#define __CELL

#include <Eigen/Core>
//#include <string>

#include "atom.hpp"

class Cell
{
private:

	Eigen::Vector3d lattice_param;
	Eigen::Vector3d lattice_angle;

	Eigen::Vector3d real_vector;
	Eigen::Vector3d reci_vector;

	int NumberOfAtoms;
	Atom* AtomList[1024];


public:

	Cell( std::string );
	
	void Show() const;

	~Cell();

};

#endif
