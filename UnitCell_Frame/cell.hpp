#ifndef __CELL
#define __CELL

#include <Eigen/Dense>	// Matrix/Vector Arithematics?
//#include <string>

#include "atom.hpp"

class Cell
{
private:

	Eigen::Vector3d lattice_param;
	Eigen::Vector3d lattice_angle;

	Eigen::Vector3d real_vector[3];		// Access? real_vector[i](0,1,2)
	Eigen::Vector3d reci_vector[3];

	double volume;

	int NumberOfAtoms;
	Atom* AtomList[1024];


public:

	Cell( std::string );
	
	void Show() const;

	virtual ~Cell();

};

#endif
