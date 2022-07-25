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

// Derivatives

Eigen::Vector3d Manager::CoulombDerivativeReal( const Cell& C, const Atom& Ai, const Atom& Aj, Eigen::Vector3d& TransVector, Eigen::Vector3d& Rij )     
{
        Eigen::Vector3d res; res.Zero();
	double r_norm = Rij.norm();
	double r_sqr  = r_norm*r_norm;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
        
        if( Ai.type == "core" && Aj.type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		res = C.TO_EV*(0.5*Ai.charge*Aj.charge)*(-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma/r_norm)-(erfc(r_norm)/C.sigma/r_sqr))/r_norm * Rij;
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
