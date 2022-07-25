#ifndef __CELL
#define __CELL

#include <Eigen/Dense>	// Matrix/Vector Arithematics?
#include <string>

#include "Atom.hpp"

#define DEF_MAX_ATOM_NUMBER 1024
#define DEF_PERIODIC_SUMM_ACCURACY 10E-25

class Cell
{

friend class Manager;

private:

	const double TO_EV = 14.39964390675221758120;

	Eigen::Vector3d lattice_param;
	Eigen::Vector3d lattice_angle;

	Eigen::Matrix3d lattice_matrix;		// lattice matrix {{a1,a2,a3},{b1,b2,b3},{c1,c2,c3}}

	Eigen::Vector3d real_vector[3];		// Accessing example : real_vector[i](0,1,2)
	Eigen::Vector3d reci_vector[3];

	double volume;				// Cell Volume

	int NumberOfAtoms;
	Atom* AtomList[DEF_MAX_ATOM_NUMBER];

	// PERIODIC SUM PARAMETERS
	double sigma;							// Sigma
	const double accuracy = DEF_PERIODIC_SUMM_ACCURACY;		// A parameter 
	double rcut,gcut;						// Cutoff tolerances
	double h_max,k_max,l_max;


	// Monopole Energy - i.e., point charge - e.g., core, shell, lone pair core ...
	double mono_real_energy, mono_reci_energy, mono_reci_self_energy;
	double mono_total_energy;


	// Cell Strain Derivatives
	Eigen::Matrix3d lattice_sd;	// ordering convention - e11(xx), e22(yy), e33(zz), e12(xy), e13(xz), e23(yz)
					//			 e1       e2       e3       e6       e5       e4

	// PERIODIC SUM WORKING VARIABLES




public:

	Cell( std::string );
	
	void ShowBasicCellInfo() const;

	void CalcCoulombEnergy();					// () field can potentially be used for adding constraints, e.g., E fields later by overloading
	void CalcCoulombDerivative();

	virtual ~Cell();

};

#endif
