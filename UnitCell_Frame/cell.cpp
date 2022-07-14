#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cstdlib>

#include "Cell.hpp"

Cell::Cell( std::string input )
{
	std::ifstream fp;

	fp.open(input);
	
	if( !fp.is_open() )
	{	std::cout << "Failed to read file ..." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		int comment_loc;
		std::string str;
		std::stringstream ss;

		for(auto str; std::get_line(fp,str); )
		{
			comment_loc = str.find('#');
			str = str.substr(0,comment_loc);
		
			if( str.length() > 0 )
			{	
				if( !str.compare("cell") )
				{


			}
	
		}
	}


}

Cell::~Cell() {}

void Cell::Show() const
{
	std::cout << lattice_vector << std::endl;
	std::cout << lattice_angle  << std::endl;
}
