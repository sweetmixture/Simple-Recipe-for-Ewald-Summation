#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cstdlib>
#include <cmath>

#include "Cell.hpp"

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
		{
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

					std::getline(fp,str);	// red one more line
					ss.clear(); ss.str("");	// clear 'ss' buffer

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

				}//	end if 'cell'

				if( !tmp.compare("atom") )
				{
					tmp.clear();
					ss >> frac[0] >> frac[1] >> frac[2] >> species >> type;
					this->AtomList[NumberOfAtoms++] = new Atom(frac[0],frac[1],frac[2],species,type);
				}

				if( !tmp.compare("species") )
				{
					tmp.clear();
					ss >> species >> type >> charge;

					for(int i=0;i<this->NumberOfAtoms;i++)
					{
						if( species == AtomList[i]->species && type == AtomList[i]->type )
						{
							AtomList[i]->charge = charge;
						}
					}
				}



			}//	end length
		}// end for getline()
	fp.close();
	}


}

Cell::~Cell()
{
	for(auto i=0;i<this->NumberOfAtoms;i++)
	{	
		delete AtomList[i];
	}
}

void Cell::Show() const
{
	using std::cout; using std::endl;
	/*
	cout << "lattice param\n";
	cout << lattice_param  << endl;
	cout << "lattice angle\n";
	cout << lattice_angle  << endl;
	*/
	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << " Basic Cell Information\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << endl;
	printf(" Cell Volume (Angs^3) : %12.6lf\n",this->volume);
	cout << endl;
	cout << "             x           y          z\n";
	printf(" a | %12.6lf%12.6lf%12.6lf\n",real_vector[0](0),real_vector[0](1),real_vector[0](2));
	printf(" b | %12.6lf%12.6lf%12.6lf\n",real_vector[1](0),real_vector[1](1),real_vector[1](2));
	printf(" c | %12.6lf%12.6lf%12.6lf\n",real_vector[2](0),real_vector[2](1),real_vector[2](2));
	cout << endl;
	cout << "             u           v          w\n";
	printf(" u | %12.6lf%12.6lf%12.6lf\n",reci_vector[0](0)/(2.*M_PI),reci_vector[0](1)/(2.*M_PI),reci_vector[0](2)/(2.*M_PI));
	printf(" v | %12.6lf%12.6lf%12.6lf\n",reci_vector[1](0)/(2.*M_PI),reci_vector[1](1)/(2.*M_PI),reci_vector[1](2)/(2.*M_PI));
	printf(" w | %12.6lf%12.6lf%12.6lf\n",reci_vector[2](0)/(2.*M_PI),reci_vector[2](1)/(2.*M_PI),reci_vector[2](2)/(2.*M_PI));
	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << " Atom fractional coordinate\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "            x           y           z       type     charge\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	for(auto i=0;i<NumberOfAtoms;i++)
	{	AtomList[i]->ShowFrac();
	}
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << " Atom Cartesian coordinate\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "            x           y           z       type     charge\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	for(auto i=0;i<NumberOfAtoms;i++)
	{	AtomList[i]->ShowFrac();
	}
}



























