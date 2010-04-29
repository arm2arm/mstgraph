#include "utils.h"
#include  <fstream>
#include <boost/lexical_cast.hpp>
using boost::lexical_cast;
using boost::bad_lexical_cast;


bool is_file_exist(const std::string &filename)
	{
	std::ifstream fin(filename.c_str(), std::ios::in);
	if(fin.is_open())
		{
		fin.close();
		return true;
		}
	fin.close();

	return false;
	}

unsigned int GetISnap(std::string str)
	{
	size_t found=str.find_last_of("_");
	unsigned int val=boost::lexical_cast<unsigned int>(str.substr(found+1));
	return val;
	}