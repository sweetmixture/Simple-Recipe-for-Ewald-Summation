#ifndef __INTERACTION_Manager
#define __INTERACTION_Manager

#include <Eigen/Dense>
#include "Cell.hpp"
#include "Atom.hpp"

class Manager // Interaction Manager
{

public:
	// Monopole - Monopole Interaction (charge vs charge)
	// decltype intrinsic - tells 'core-core', 'shell-core', 'core-shgl'
	double CoulombMonoMonoReal( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij );
	double CoulombMonoMonoSelf( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij );
	double CoulombMonoMonoReci( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij );

	// First Geometric Derivatives
	Eigen::Vector3d CoulombDerivativeReal( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij );
	Eigen::Vector3d CoulombDerivativeReci( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij );

	// Strain Derivatives
	Eigen::Matrix3d StrainDerivativeReal( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij );
	Eigen::Matrix3d StrainDerivativeReci( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij );


};

#endif
