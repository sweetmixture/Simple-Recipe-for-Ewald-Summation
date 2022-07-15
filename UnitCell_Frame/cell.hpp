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

public:

	Cell( std::string );
	
	void Show() const;

	~Cell();

};

#endif
