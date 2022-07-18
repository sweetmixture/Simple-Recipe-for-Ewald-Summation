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
		printf("%12.4lf%12.4lf%12.4lf\n",frac_x,frac_y,frac_z);
		frac(0) = frac_x;
		frac(1) = frac_y;
		frac(2) = frac_z;
		std::cout << "Eigen out\n";
		std::cout << frac << std::endl;		

		this->species = species;
		this->type    = type;
		std::cout << "Atom Init\n";
	}

	virtual void Show() const
	{
		std::cout << frac << std::endl;
		std::cout << "## Atom Frac\n";
		std::cout << "frac : " << frac(0) << " " << frac(1) << " " << frac(2) << " " << species << " " << type << std::endl;
	}

	virtual ~Atom()								// Explicit Destructor Check
	{	
		cnt++;
		std::cout << "Destructor ~" << cnt << std::endl;
	}

};

#endif
