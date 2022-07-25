#include "Manager.hpp"

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

double Manager::CoulombMonoMonoReal( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij )
{
        double res;
        // TransVector = h*a + k*b + l*c
        // Rij         = Ai.r - Aj.r - TransVector;
        
        if( Ai.type == "core" && Aj.type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {       res = 0.5*(Ai.charge*Aj.charge)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
                return res;
        }       
        
        if( Ai.type == "core" && Aj.type == "shel" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "core" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "shel" )
        {
                return res;
        }       
        
        return 0;
}       

double Manager::CoulombMonoMonoSelf( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij )
{
        double res;
        // TransVector = h*a + k*b + l*c
        // Rij         = Ai.r - Aj.r - TransVector;
        
        if( Ai.type == "core" && Aj.type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {       res = -(Ai.charge*Aj.charge)/C.sigma/sqrt(M_PI) * C.TO_EV;
                return res;
        }       
        
        if( Ai.type == "core" && Aj.type == "shel" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "core" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "shel" )
        {
                return res;
        }       
        
        return 0;
}       

double Manager::CoulombMonoMonoReci( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij )
{
        double res;
	double g_norm = TransVector.norm();
	double g_sqr  = g_norm*g_norm;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
        
        if( Ai.type == "core" && Aj.type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {	res = (2.*M_PI/C.volume)*(Ai.charge*Aj.charge) \
			* exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
                return res;
        }       
        
        if( Ai.type == "core" && Aj.type == "shel" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "core" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "shel" )
        {
                return res;
        }       
        
        return 0;
} 

// Geometric (RAW) Derivatives

Eigen::Vector3d Manager::CoulombDerivativeReal( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij )     
{
        Eigen::Vector3d res; res.Zero();
	double r_norm = Rij.norm();
	double r_sqr  = r_norm*r_norm;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
        
        if( Ai.type == "core" && Aj.type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		res = C.TO_EV*(0.5*Ai.charge*Aj.charge)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm * Rij;
                return res;
        }       
        
        if( Ai.type == "core" && Aj.type == "shel" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "core" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "shel" )
        {
                return res;
        }       
        
        return res;
}       

Eigen::Vector3d Manager::CoulombDerivativeReci( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij )     
{
        Eigen::Vector3d res; res.Zero();
	double g_norm = TransVector.norm();
	double g_sqr  = g_norm*g_norm;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
        
        if( Ai.type == "core" && Aj.type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		res = C.TO_EV*((2.*M_PI)/C.volume)*Ai.charge*Aj.charge*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij) * TransVector;
                return res;
        }       
        
        if( Ai.type == "core" && Aj.type == "shel" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "core" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "shel" )
        {
                return res;
        }       
        
        return res;
}       

// Strain Derivative

Eigen::Matrix3d Manager::StrainDerivativeReal( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij )
{
        Eigen::Matrix3d res; res.Zero();
	double r_norm = Rij.norm();
	double r_sqr  = r_norm*r_norm;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
        
        if( Ai.type == "core" && Aj.type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		res(0,0) = C.TO_EV*(0.5*Ai.charge*Aj.charge)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr)) * Rij(0)/r_norm * Rij(0);
		res(1,1) = C.TO_EV*(0.5*Ai.charge*Aj.charge)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr)) * Rij(1)/r_norm * Rij(1);
		res(2,2) = C.TO_EV*(0.5*Ai.charge*Aj.charge)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr)) * Rij(2)/r_norm * Rij(2);
		res(1,2) = C.TO_EV*(0.5*Ai.charge*Aj.charge)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr)) * Rij(1)/r_norm * Rij(2);
		res(0,2) = C.TO_EV*(0.5*Ai.charge*Aj.charge)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr)) * Rij(0)/r_norm * Rij(2);
		res(0,1) = C.TO_EV*(0.5*Ai.charge*Aj.charge)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr)) * Rij(0)/r_norm * Rij(1);

                return res;
        }       
        
        if( Ai.type == "core" && Aj.type == "shel" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "core" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "shel" )
        {
                return res;
        }       
        
        return res;
}       

Eigen::Matrix3d Manager::StrainDerivativeReci( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij )
{
        Eigen::Matrix3d res;
	double g_norm = TransVector.norm();
	double g_sqr  = g_norm*g_norm;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
        
        if( Ai.type == "core" && Aj.type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
		res(0,0) = C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij) * TransVector(0) * Rij(0);
		res(1,1) = C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij) * TransVector(1) * Rij(1);
		res(2,2) = C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij) * TransVector(2) * Rij(2);
		res(1,2) = C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij) * TransVector(1) * Rij(2);
		res(0,2) = C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij) * TransVector(0) * Rij(2);
		res(0,1) = C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij) * TransVector(0) * Rij(1);

		// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
		res(0,0) += C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*(-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(0)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)
			-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(0)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij)-exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*Rij(0)*sin(TransVector.adjoint()*Rij))*-TransVector(0);
		res(1,1) += C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*(-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(1)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)
			-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(1)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij)-exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*Rij(1)*sin(TransVector.adjoint()*Rij))*-TransVector(1);
		res(2,2) += C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*(-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(2)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)
			-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(2)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij)-exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*Rij(2)*sin(TransVector.adjoint()*Rij))*-TransVector(2);
		res(1,2) += C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*(-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(2)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)
			-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(2)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij)-exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*Rij(2)*sin(TransVector.adjoint()*Rij))*-TransVector(1);
		res(0,2) += C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*(-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(2)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)
			-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(2)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij)-exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*Rij(2)*sin(TransVector.adjoint()*Rij))*-TransVector(0);
		res(0,1) += C.TO_EV*((2.*M_PI)/C.volume)*(Ai.charge*Aj.charge)*(-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(1)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)
			-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)*TransVector(1)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij)-exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*Rij(1)*sin(TransVector.adjoint()*Rij))*-TransVector(0);

		// Strain derivative (3) - w.r.t cell volume in the reciprocal space
		res(0,0) += -C.TO_EV*((2.*M_PI)/C.volume/C.volume)*C.volume*(Ai.charge*Aj.charge)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);
		res(1,1) += -C.TO_EV*((2.*M_PI)/C.volume/C.volume)*C.volume*(Ai.charge*Aj.charge)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);
		res(2,2) += -C.TO_EV*((2.*M_PI)/C.volume/C.volume)*C.volume*(Ai.charge*Aj.charge)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);

                return res;
        }       
        
        if( Ai.type == "core" && Aj.type == "shel" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "core" )
        {
                return res;
        }       
        
        if( Ai.type == "shel" && Aj.type == "shel" )
        {
                return res;
        }       
        
        return res;
}       
