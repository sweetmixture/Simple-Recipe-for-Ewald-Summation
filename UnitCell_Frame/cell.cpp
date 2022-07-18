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
		std::string str, tmp;
		std::stringstream ss;

		double frac[3];
		std::string species;
		std::string type;


		for(std::string str; std::getline(fp,str); )
		{
			comment_loc = str.find('#');
			str = str.substr(0,comment_loc);
		
			if( str.length() > 0 )
			{
				if( !str.compare("cell") )
				{
					std::getline(fp,str);
					ss << str;
					
					ss >> this->lattice_param(0);
					ss >> this->lattice_param(1);
					ss >> this->lattice_param(2);

					ss >> this->lattice_angle(0);
					ss >> this->lattice_angle(1);
					ss >> this->lattice_angle(2);
				}
				
				if( !str.compare("fractional") )
				{
					for(std::string str; std::getline(fp,str); )
					{
						comment_loc = str.find('#');
						str = str.substr(0,comment_loc);
						ss.clear();	ss.str("");
						ss << str;
						ss >> tmp;

						//std::cout << str << std::endl;
						//std::cout << str.length() << std::endl;
					

						if( !tmp.compare("atom") )
						{
							tmp.clear();	
							ss >> frac[0] >> frac[1] >> frac[2] >> species >> type;

							printf("%12.6lf%12.6lf%12.6lf%12.6s%12.6s\n", frac[0],frac[1],frac[2],species.c_str(),type.c_str());
							this->AtomList[NumberOfAtoms++] = new Atom(frac[0],frac[1],frac[2],species,type);
						}
					}
					std::cout << NumberOfAtoms << std::endl;
					break;
				}
			}
		}
		
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

	cout << "lattice param\n";
	cout << lattice_param  << endl;
	cout << "lattice angle\n";
	cout << lattice_angle  << endl;


	for(auto i=0;i<NumberOfAtoms;i++)
	{
		AtomList[i]->Show();
	}
}
