#ifndef __ATOM
#define __ATOM

#include <Eigen/Core>
#include <iostream>
#include <string>

static int cnt = 0;

class Atom
{

friend class Cell;				// Allowing access previlege to 'Cell' class, i.e., Cell can access privates of Atom
friend class Manager;

private:

	double charge;				// Atom charge

	std::string species;			// Atom name
	std::string type;			// Atom description - core, shel, lone (lonepair)

	Eigen::Vector3d frac;			// fractional coord
	Eigen::Vector3d cart;			// fractional -> Cartesian

	Eigen::Vector3d cart_gd;		// g.d. -> geometric derivative
	Eigen::Vector3d cart_gd_int;		// g.d. -> internal ; same as what GULP gives

public:

	Atom( const double frac_x, const double frac_y, const double frac_z, std::string species, std::string type, const Eigen::Vector3d* lattice_vector )
	{
		frac(0) = frac_x;
		frac(1) = frac_y;
		frac(2) = frac_z;

		this->species = species;
		this->type    = type;

		cart = frac(0)*lattice_vector[0] + frac(1)*lattice_vector[1] + frac(2)*lattice_vector[2];
		// cartesian r = x_frac * a + y_frac * b + z_frac * c;
	}

	virtual void SetFeature( const double charge )
	{
		this->charge = charge;
	}

	virtual void ShowFrac() const
	{
		printf("%4.3s%12.6lf%12.6lf%12.6lf%8.4s%12.6lf\n",species.c_str(),frac(0),frac(1),frac(2),type.c_str(),charge);
	}

	virtual void ShowCart() const
	{
		printf("%4.3s%12.6lf%12.6lf%12.6lf%8.4s%12.6lf\n",species.c_str(),cart(0),cart(1),cart(2),type.c_str(),charge);
	}

	virtual ~Atom()								// Explicit Destructor Check
	{	
		cnt++;
		std::cout << "Atom Destructor ~" << cnt << std::endl;
	}

};

class Shell : public Atom
{
private:

	double shell_charge;			// Shell charge

	Eigen::Vector3d shell_frac;		// Below all same but for 'Shell'
	Eigen::Vector3d shell_cart;
	
	Eigen::Vector3d shell_cart_gd;
	Eigen::Vector3d shell_cart_gd_int;

public:

};

class LonePair : public Atom
{



};

#endif
