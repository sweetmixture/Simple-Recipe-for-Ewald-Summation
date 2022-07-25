#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include <cstdlib>
#include <cmath>

#include "cell.hpp"
#include "Manager.hpp"

static int line_cnt = 0;

Cell::Cell( std::string input )
{
	this->NumberOfAtoms = 0;

	std::ifstream fp;

	fp.open(input);
	
	if( !fp.is_open() )
	{	std::cout << "Failed to read file ..." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		int comment_loc;		// Comment locations 
		std::string str, tmp;		//
		std::stringstream ss;

		double frac[3];			// tmp load factional coordinate of atoms
		std::string species;		// species - e.g., Mg, O etc ...
		std::string type;		// type    - shell, core ?.. (type of describing a species)
		double charge;


		double a,b,c,al,be,ga;		// worksapce - real/reci lattice vector
		double v[3];			// workspace - real/reci lattice vector
	
		Eigen::Vector3d vv;		// workspace - vector operation

		for(std::string str; std::getline(fp,str); )
		{	line_cnt++;

			comment_loc = str.find('#');
			str = str.substr(0,comment_loc);
			tmp.clear();

			if( str.length() > 0 )
			{
				ss.clear(); ss.str("");
				ss << str;
				ss >> tmp;		

				if( !tmp.compare("cell") )
				{
					tmp.clear();

					std::getline(fp,str);	line_cnt++; 	// red one more line
					ss.clear(); ss.str("");			// clear 'ss' buffer

					ss << str;
					
					ss >> this->lattice_param(0);
					ss >> this->lattice_param(1);
					ss >> this->lattice_param(2);

					ss >> this->lattice_angle(0);	this->lattice_angle(0) = this->lattice_angle(0)/180.*M_PI;
					ss >> this->lattice_angle(1);	this->lattice_angle(1) = this->lattice_angle(1)/180.*M_PI;
					ss >> this->lattice_angle(2);	this->lattice_angle(2) = this->lattice_angle(2)/180.*M_PI;	// angle 'degree' -> 'radian'

					// PROCESSING REAL LATTICE VECTORS / http://gisaxs.com/index.php/Unit_cell
					v[0] = this->lattice_param(0);
					v[1] = 0.;
					v[2] = 0.;
					this->real_vector[0] << v[0],v[1],v[2];			// set real lattice vector a = (L1,0,0);

					v[0] = this->lattice_param(1)*std::cos(this->lattice_angle(2));
					v[1] = this->lattice_param(1)*std::sin(this->lattice_angle(2));
					v[2] = 0.;
					this->real_vector[1] << v[0],v[1],v[2];			// set real lattice vector b = (L2*Cos(gamma),L2*Sin(gamma),0);

					v[0] = this->lattice_param(2)*std::cos(this->lattice_angle(1));
					v[1] = this->lattice_param(2)*(std::cos(this->lattice_angle(0))-std::cos(this->lattice_angle(1))*std::cos(this->lattice_angle(2)))/std::sin(this->lattice_angle(2));
					v[2] = this->lattice_param(2)*std::sqrt(1.-pow(std::cos(this->lattice_angle(1)),2.)-pow(v[1]/this->lattice_param(2),2.));		
					this->real_vector[2] << v[0],v[1],v[2];			// set real lattice vector c

					// Setting LatticeMatrix - Eigen::Matrix3d

					this->lattice_matrix << this->real_vector[0](0),this->real_vector[0](1),this->real_vector[0](2), // order	{a1,a2,a3}
								this->real_vector[1](0),this->real_vector[1](1),this->real_vector[1](2), //		{b1,b2,b3}
								this->real_vector[2](0),this->real_vector[2](1),this->real_vector[2](2); //		{c1,c2,c3}

					// Expected result : T_hkl = a h + b k + c l		
					//
					//

					// PROCESSING RECIPROCAL LATTICE VECTORS
					vv = this->real_vector[1].cross(this->real_vector[2]);	// get b x c
					this->volume = this->real_vector[0].adjoint()*vv;	// get a.(b x c)

					// u vector
					vv = 2.*M_PI*this->real_vector[1].cross(this->real_vector[2])/this->volume;	// (b x c) / Vcell
					this->reci_vector[0] << vv(0), vv(1), vv(2);
					// v vector
					vv = 2.*M_PI*this->real_vector[2].cross(this->real_vector[0])/this->volume;	// (c x a) / Vcell
					this->reci_vector[1] << vv(0), vv(1), vv(2);
					// w vector
					vv = 2.*M_PI*this->real_vector[0].cross(this->real_vector[1])/this->volume;	// (a x b) / Vcell
					this->reci_vector[2] << vv(0), vv(1), vv(2);

					if( this->volume < 0. ) { this->volume = this->volume*-1.; }

					// Expected result : G_hkl = 2PI u h + 2PI v k + 2PI w l	// Note that 2.*M_PI factor is multiplied
					//
					//

				}// end if 'cell'

				if( !tmp.compare("atom") ) // DynamicMemory Mangement - Atom
				{
					tmp.clear();
					ss >> frac[0] >> frac[1] >> frac[2] >> species >> type;

					if( !type.compare("core") )
					{	this->AtomList[NumberOfAtoms++] = new Atom(frac[0],frac[1],frac[2],species,type,real_vector);	// here 'real_vector' contains lattice vectors
					}
					else if( !type.compare("shel") )
					{	// Need appropriate initialiser
					}
					else if( !type.compare("lone") )
					{	// Need appropriate initialiser
					}
					else	// TYPE CHECK
					{	std::cout << "error found in atom type ... ";
						std::cout << "specified type : " << "'" << type << "'" << " is not found in " << "'" << input << "'\n"; 
						std::cout << "list of allowed types : 'core', 'shel', and 'lone'\n";
						std::cout << "<SEE LINE> : " << line_cnt << std::endl;
						std::exit(EXIT_FAILURE);
					}
				}// end if 'atom'

				if( !tmp.compare("species") )
				{
					tmp.clear();
					ss >> species >> type >> charge;

					for(int i=0;i<this->NumberOfAtoms;i++)
					{	
						if( species == AtomList[i]->species && type == AtomList[i]->type )
						{
							if( !type.compare("core") )
							{	AtomList[i]->SetFeature(charge);
								/*
								double k;
								ss >> k;
								std::cout << k << "test\n";
								*/	// if 'stringstream' ss is empty, gives noting ... expected output '0'
							}
							else if( !type.compare("shel") )
							{	// Need appropriate initialiser
							}
							else if( !type.compare("lone") )
							{	// Need appropriate initialiser
							}
						}
					}

					// TYPE_CEHCK .. need implementation
					/*
						else	// TYPE CHECK
						{	std::cout << "no matching types found in species with atoms ... ";
							std::cout << "specified species " << "'" << species << "' ";
							std::cout << "or type " << "'" << type << "'" << " is not found in " << "'" << input << "' atom list\n"; 
							std::cout << "<SEE LINE> : " << line_cnt << std::endl;
							std::cout << "list of allowed types : 'core', 'shel', and 'lone'\n";
							std::exit(EXIT_FAILURE);
						}
					*/
				}// end if 'species'






			}//	end length
		}// end for getline()
	fp.close();
	}

	// Prepare Parameters - Periodic Summation
	this->sigma = pow(NumberOfAtoms*M_PI*M_PI*M_PI/this->volume/this->volume,-1./6.);
	this->rcut  = std::sqrt(-log(this->accuracy)*this->sigma*this->sigma);
	this->gcut  = 2./this->sigma*std::sqrt(-log(this->accuracy));
	
	auto lattice_min = std::min({this->real_vector[0].norm(),this->real_vector[1].norm(),this->real_vector[2].norm(),
				     this->reci_vector[0].norm(),this->reci_vector[1].norm(),this->reci_vector[2].norm()});
	int max_grid     = static_cast<int>(std::max(this->rcut,this->gcut)/lattice_min+1);
	this->h_max = this->k_max = this->l_max = max_grid; 

//std::cout << "lattice min : " << lattice_min << "\n";
//std::cout << "max grid : " << max_grid << std::endl;
//std::cout << "rcut/gcut: " << rcut << "/" << gcut << std::endl;
//std::cout << "sigma    : " << sigma << std::endl;

}	// END 'Cell' class initilising


void Cell::CalcCoulombEnergy()
{
using std::cos, std::sin, std::sqrt;

	Manager manager;	// Managing class - interaction

	Eigen::Vector3d trans;
	double trans_norm;

	Eigen::Vector3d delta_r;
	Eigen::Vector3d delta_g;
	double g_norm, g_sqr;
	double r_norm, r_sqr;		// Aux variables

	this->mono_real_energy = this->mono_reci_energy = this->mono_reci_self_energy = this->mono_total_energy = 0.;	// Initialising energies

	for(int i=0;i<this->NumberOfAtoms;i++)
	{
		for(int j=0;j<this->NumberOfAtoms;j++)
		{
			for(int h = -this->h_max ; h <= this->h_max ; h++)
			{
				for(int k = -this->k_max ; k <= this->k_max ; k++)
				{
					for(int l = -this->l_max ; l <= this->l_max ; l++)
					{
						//// START REAL SPACE
						trans   = h*this->real_vector[0] + k*this->real_vector[1] + l*this->real_vector[2];	// T = h*a + k*b +l*c
						delta_r = this->AtomList[i]->cart - this->AtomList[j]->cart - trans;			// ri - rj - T
						r_norm  = delta_r.norm();
						r_sqr   = r_norm*r_norm;

						if( (trans.norm()) < this->rcut )
						{	
							if( h == 0 && k == 0 && l == 0 )
							{	
							    if( i != j )
							    {
								this->mono_real_energy += manager.CoulombMonoMonoReal(*this,*this->AtomList[i],*this->AtomList[j],trans,delta_r);
							    }	// h=k=l=0 (central image) - excluding self interaction
							}
							else
							{
								this->mono_real_energy += manager.CoulombMonoMonoReal(*this,*this->AtomList[i],*this->AtomList[j],trans,delta_r);
							}
						}
						//// END REAL SPACE

						//## START RECIPROCAL SPACE
						trans   = h*this->reci_vector[0] + k*this->reci_vector[1] + l*this->reci_vector[2];	// G = h*2pi*u + k*2pi*v + l*2pi*w
						g_norm  = trans.norm();									// |G|
						g_sqr   = g_norm*g_norm;								// |G|^2

						delta_r = this->AtomList[i]->cart - this->AtomList[j]->cart;				// ri - rj
						
						if( (trans.norm()) < this->gcut )
						{
							if( h == 0 && k == 0 && l == 0 )
							{
							    if( i == j )	// Self interaction
							    {	
								this->mono_reci_self_energy += manager.CoulombMonoMonoSelf(*this,*this->AtomList[i],*this->AtomList[j],trans,delta_r);
							    }
							}
							else
							{	
								this->mono_reci_self_energy += manager.CoulombMonoMonoSelf(*this,*this->AtomList[i],*this->AtomList[j],trans,delta_r);
							}
						}
					}// end l
				}// end k
			}//end h
		}// end j
	}// end i

	this->mono_total_energy = this->mono_real_energy + this->mono_reci_energy + this->mono_reci_self_energy;


std::cout << "\n";
std::cout << "Output Val\n";
std::cout << "Real : " << std::setprecision(16) << this->mono_real_energy << std::endl;
std::cout << "Reci : " << std::setprecision(16) << this->mono_reci_energy << std::endl;
std::cout << "Self : " << std::setprecision(16) << this->mono_reci_self_energy << std::endl;
std::cout << "\n";
std::cout << "Real : " << std::setprecision(16) << this->mono_real_energy << std::endl;
std::cout << "Rc+S : " << std::setprecision(16) << this->mono_reci_self_energy + this->mono_reci_energy << std::endl;
std::cout << "Tota : " << std::setprecision(16) << this->mono_total_energy << std::endl;
std::cout << "\n";

	return;
}

void Cell::CalcCoulombDerivative()
{

/*
	Calculates Derivatives

	1) Raw Geometric
	2) Internal Geometric
	3) Strain
*/
	Manager manager;

	Eigen::Vector3d trans;
	double trans_norm;

	Eigen::Vector3d delta_r;
	Eigen::Vector3d delta_g;
	double g_norm, g_sqr;
	double r_norm, r_sqr;		// Aux variables

	double q_prod;

	// Init derivative fields
	for(int i=0;i<this->NumberOfAtoms;i++)
	{
		this->AtomList[i]->cart_gd.Zero();
	}	
	this->lattice_sd.Zero();	// lattice strain derivatives

	for(int i=0;i<this->NumberOfAtoms;i++)
	{
		for(int j=0;j<this->NumberOfAtoms;j++)
		{
			q_prod = this->AtomList[i]->charge*this->AtomList[j]->charge;			

			for(int h = -this->h_max ; h <= this->h_max ; h++)
			{
				for(int k = -this->k_max ; k <= this->k_max ; k++)
				{
					for(int l = -this->l_max ; l <= this->l_max ; l++)
					{

						//// START REAL SPACE
						trans   = h*this->real_vector[0] + k*this->real_vector[1] + l*this->real_vector[2];	// T = h*a + k*b +l*c
						delta_r = this->AtomList[i]->cart - this->AtomList[j]->cart - trans;			// ri - rj - T
						r_norm  = delta_r.norm();
						r_sqr   = r_norm*r_norm;

						if( (trans.norm()) < this->rcut )
						{	
							if( h == 0 && k == 0 && l == 0 )
							{	
							    if( i != j )
							    {	// Raw Geometric Derivatives
								this->AtomList[i]->cart_gd += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr))/r_norm * delta_r;
								this->AtomList[j]->cart_gd -= this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr))/r_norm * delta_r;
								//this->AtomList[i]->cart_gd += manager.CoulombDerivativeReal(*this,*this->AtomList[i],*this->AtomList[j],trans,delta_r);
								//this->AtomList[j]->cart_gd -= manager.CoulombDerivativeReal(*this,*this->AtomList[i],*this->AtomList[j],trans,delta_r);

								// Strain Derivatives
								this->lattice_sd(0,0) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(0)/r_norm * delta_r(0);
								this->lattice_sd(1,1) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(1)/r_norm * delta_r(1);
								this->lattice_sd(2,2) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(2)/r_norm * delta_r(2);
								this->lattice_sd(1,2) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(1)/r_norm * delta_r(2);
								this->lattice_sd(0,2) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(0)/r_norm * delta_r(2);
								this->lattice_sd(0,1) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(0)/r_norm * delta_r(1);

							    }	// h=k=l=0 (central image) - excluding self interaction
							}
							else
							{	// Raw Geometric Derivatives
								this->AtomList[i]->cart_gd += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr))/r_norm * delta_r;
								this->AtomList[j]->cart_gd -= this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr))/r_norm * delta_r;

								// Strain Derivatives
								this->lattice_sd(0,0) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(0)/r_norm * delta_r(0);
								this->lattice_sd(1,1) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(1)/r_norm * delta_r(1);
								this->lattice_sd(2,2) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(2)/r_norm * delta_r(2);
								this->lattice_sd(1,2) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(1)/r_norm * delta_r(2);
								this->lattice_sd(0,2) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(0)/r_norm * delta_r(2);
								this->lattice_sd(0,1) += this->TO_EV*(0.5*q_prod)*((-2./this->sigma/sqrt(M_PI))*(exp(-r_sqr/this->sigma/this->sigma)/r_norm)-(erfc(r_norm/this->sigma)/r_sqr)) * delta_r(0)/r_norm * delta_r(1);
							}
						}
						//// END REAL SPACE

						//## START RECIPROCAL SPACE
						trans   = h*this->reci_vector[0] + k*this->reci_vector[1] + l*this->reci_vector[2];	// G = h*2pi*u + k*2pi*v + l*2pi*w
						g_norm  = trans.norm();									// |G|
						g_sqr   = g_norm*g_norm;								// |G|^2

						delta_r = this->AtomList[i]->cart - this->AtomList[j]->cart;				// ri - rj
						
						if( (trans.norm()) < this->gcut )
						{
							if( h == 0 && k == 0 && l == 0 )
							{
								if( i == j )	// Self interaction
								{
									// Do Nothing ... since the derivatives are zero
								}
							}
							else
							{	// Raw Geometric Derivative
								this->AtomList[i]->cart_gd += this->TO_EV*((2.*M_PI)/this->volume)*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * -sin(trans.adjoint()*delta_r) * trans;
								this->AtomList[j]->cart_gd -= this->TO_EV*((2.*M_PI)/this->volume)*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * -sin(trans.adjoint()*delta_r) * trans;

								// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
								this->lattice_sd(0,0) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * -sin(trans.adjoint()*delta_r) * trans(0) * delta_r(0);
								this->lattice_sd(1,1) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * -sin(trans.adjoint()*delta_r) * trans(1) * delta_r(1);
								this->lattice_sd(2,2) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * -sin(trans.adjoint()*delta_r) * trans(2) * delta_r(2);
								this->lattice_sd(1,2) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * -sin(trans.adjoint()*delta_r) * trans(1) * delta_r(2);
								this->lattice_sd(0,2) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * -sin(trans.adjoint()*delta_r) * trans(0) * delta_r(2);
								this->lattice_sd(0,1) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * -sin(trans.adjoint()*delta_r) * trans(0) * delta_r(1);

								// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
								this->lattice_sd(0,0) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod
										* (-2.*exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(0)/g_sqr/g_sqr*cos(trans.adjoint()*delta_r)
										- 0.5 *exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(0)/g_sqr*this->sigma*this->sigma*cos(trans.adjoint()*delta_r)
										-      exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr*delta_r(0)*sin(trans.adjoint()*delta_r)) * -trans(0);
								this->lattice_sd(1,1) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod
										* (-2.*exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(1)/g_sqr/g_sqr*cos(trans.adjoint()*delta_r)
										- 0.5 *exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(1)/g_sqr*this->sigma*this->sigma*cos(trans.adjoint()*delta_r)
										-      exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr*delta_r(1)*sin(trans.adjoint()*delta_r)) * -trans(1);
								this->lattice_sd(2,2) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod
										* (-2.*exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(2)/g_sqr/g_sqr*cos(trans.adjoint()*delta_r)
										- 0.5 *exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(2)/g_sqr*this->sigma*this->sigma*cos(trans.adjoint()*delta_r)
										-      exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr*delta_r(2)*sin(trans.adjoint()*delta_r)) * -trans(2);
								this->lattice_sd(1,2) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod
										* (-2.*exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(2)/g_sqr/g_sqr*cos(trans.adjoint()*delta_r)
										- 0.5 *exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(2)/g_sqr*this->sigma*this->sigma*cos(trans.adjoint()*delta_r)
										-      exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr*delta_r(2)*sin(trans.adjoint()*delta_r)) * -trans(1);
								this->lattice_sd(0,2) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod
										* (-2.*exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(2)/g_sqr/g_sqr*cos(trans.adjoint()*delta_r)
										- 0.5 *exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(2)/g_sqr*this->sigma*this->sigma*cos(trans.adjoint()*delta_r)
										-      exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr*delta_r(2)*sin(trans.adjoint()*delta_r)) * -trans(0);
								this->lattice_sd(0,1) += this->TO_EV*((2.*M_PI)/this->volume)*q_prod
										* (-2.*exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(1)/g_sqr/g_sqr*cos(trans.adjoint()*delta_r)
										- 0.5 *exp(-0.25*this->sigma*this->sigma*g_sqr)*trans(1)/g_sqr*this->sigma*this->sigma*cos(trans.adjoint()*delta_r)
										-      exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr*delta_r(1)*sin(trans.adjoint()*delta_r)) * -trans(0);

								// Strain derivative (3) - w.r.t cell volume in the reciprocal space
								this->lattice_sd(0,0) += -this->TO_EV*((2.*M_PI)/this->volume/this->volume)*this->volume*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * cos(trans.adjoint()*delta_r);
								this->lattice_sd(1,1) += -this->TO_EV*((2.*M_PI)/this->volume/this->volume)*this->volume*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * cos(trans.adjoint()*delta_r);
								this->lattice_sd(2,2) += -this->TO_EV*((2.*M_PI)/this->volume/this->volume)*this->volume*q_prod*exp(-0.25*this->sigma*this->sigma*g_sqr)/g_sqr * cos(trans.adjoint()*delta_r);
							}
						}
					}// end l
				}// end k
			}//end h
		}// end j
	}// end i


	// Convert raw geometric derivative -> internal geometric derivative
	for(int i=0;i<this->NumberOfAtoms;i++)
	{	AtomList[i]->cart_gd_int = this->lattice_matrix * AtomList[i]->cart_gd;
	}
	// Symmeterise lattice strain derivative matrix
	this->lattice_sd(1,0) = this->lattice_sd(0,1);
	this->lattice_sd(2,0) = this->lattice_sd(0,2);
	this->lattice_sd(2,1) = this->lattice_sd(1,2);

	std::cout << "1> Law Cartesian Derivative Check\n";
	for(int i=0;i<this->NumberOfAtoms;i++)
	{	
		printf("%12.4s\t%12.6lf\t%12.6lf\t%12.6lf\n",this->AtomList[i]->species.c_str(),this->AtomList[i]->cart_gd(0),this->AtomList[i]->cart_gd(1),this->AtomList[i]->cart_gd(2));
	}
	std::cout << "2> Internal Cartesian Derivative Check\n";
	for(int i=0;i<this->NumberOfAtoms;i++)
	{	
		printf("%12.4s\t%12.6lf\t%12.6lf\t%12.6lf\n",this->AtomList[i]->species.c_str(),this->AtomList[i]->cart_gd_int(0),this->AtomList[i]->cart_gd_int(1),this->AtomList[i]->cart_gd_int(2));
	}
	std::cout << "3> Strain Derivative Check\n";
	for(int i=0;i<3;i++)
	{	for(int j=0;j<3;j++)
		{
			printf("%12.6lf\t",this->lattice_sd(i,j));
		}
		std::cout << "\n";
	}
	
	return;
}

Cell::~Cell()
{
	for(auto i=0;i<this->NumberOfAtoms;i++)
	{	
		delete AtomList[i];
	}
}



























////	////	////	////

////	 OUTPUT CHECK

////	////	////	////

void Cell::ShowBasicCellInfo() const
{

using std::cout;
using std::endl;

	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << " Basic Cell Information\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << endl;
	printf("     Cell Volume (Angs^3) : %12.6lf\n",this->volume);
	cout << endl;
	cout << "                    x           y          z\n";
	printf("        a | %12.6lf%12.6lf%12.6lf\n",real_vector[0](0),real_vector[0](1),real_vector[0](2));
	printf("        b | %12.6lf%12.6lf%12.6lf\n",real_vector[1](0),real_vector[1](1),real_vector[1](2));
	printf("        c | %12.6lf%12.6lf%12.6lf\n",real_vector[2](0),real_vector[2](1),real_vector[2](2));
	cout << endl;
	cout << "                    u           v          w\n";
	printf("        u | %12.6lf%12.6lf%12.6lf\n",reci_vector[0](0)/(2.*M_PI),reci_vector[0](1)/(2.*M_PI),reci_vector[0](2)/(2.*M_PI));
	printf("        v | %12.6lf%12.6lf%12.6lf\n",reci_vector[1](0)/(2.*M_PI),reci_vector[1](1)/(2.*M_PI),reci_vector[1](2)/(2.*M_PI));
	printf("        w | %12.6lf%12.6lf%12.6lf\n",reci_vector[2](0)/(2.*M_PI),reci_vector[2](1)/(2.*M_PI),reci_vector[2](2)/(2.*M_PI));
	cout << endl;
	printf("  2pi * u | %12.6lf%12.6lf%12.6lf\n",reci_vector[0](0),reci_vector[0](1),reci_vector[0](2));
	printf("  2pi * v | %12.6lf%12.6lf%12.6lf\n",reci_vector[1](0),reci_vector[1](1),reci_vector[1](2));
	printf("  2pi * w | %12.6lf%12.6lf%12.6lf\n",reci_vector[2](0),reci_vector[2](1),reci_vector[2](2));
	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "    Atom fractional coordinate\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "            x           y           z       type     charge\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	for(auto i=0;i<NumberOfAtoms;i++)
	{	AtomList[i]->ShowFrac();
	}
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "    Atom Cartesian coordinate\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "            x           y           z       type     charge\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	for(auto i=0;i<NumberOfAtoms;i++)
	{	AtomList[i]->ShowCart();
	}
}
