#ifndef __ATOM
#define __ATOM

#include <Eigen/Core>
#include <string>

static int cnt = 0;

class Atom
{

friend class Cell;				// Allowing access to the Cell class - unit cell information

private:

	double charge;				// Atom charge

	std::string species;			// Atom name
	std::string type;

	Eigen::Vector3d frac;			// fractional coord
	Eigen::Vector3d cart;			// fractional -> Cartesian

	Eigen::Vector3d cart_gd;		// g.d. -> geometric derivative
	Eigen::Vector3d cart_gd_int;		// g.d. -> internal ; same as what GULP gives

public:

	Atom( const double frac_x, const double frac_y, const double frac_z, std::string species, std::string type )
	{
		frac(0) = frac_x;
		frac(1) = frac_y;
		frac(2) = frac_z;

		this->species = species;
		this->type    = type;
	}

	virtual void ShowFrac() const
	{
		printf("%4.3s%12.6lf%12.6lf%12.6lf%8.4s%12.6lf\n",species.c_str(),frac(0),frac(1),frac(2),type.c_str(),charge);
	}

	virtual void ShowReal() const
	{

	}

	virtual ~Atom()								// Explicit Destructor Check
	{	
		cnt++;
		std::cout << "Destructor ~" << cnt << std::endl;
	}

};

#endif
