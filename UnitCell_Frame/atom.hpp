#ifndef __ATOM
#define __ATOM

#include <Eigen/Core>
#include <string>

class Atom
{

private:

	std::string tag;			// Atom name
	double q;				// Atom charge

	Eigen::Vector3d frac;			// fractional coord
	Eigen::Vector3d cart;			// fractional -> Cartesian

	Eigen::Vector3d cart_force;		// cumulative force

public:

};

#endif
