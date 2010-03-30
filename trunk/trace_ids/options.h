#ifndef _MYOPTIONS_
#define _MYOPTIONS_
#include <iostream>
#include <exception>
#include <fstream>
#include <cstdlib>
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ostream_iterator;

class COptions{
public:
	COptions(int argc, char* argv[]);
	const int is_bad(){return m_status;};
//////////////////////////////////////////
	vector<string> m_snapshotList;
	vector<int> m_IDlist;
	string m_file_out;
	int m_type;
private:
	int m_status;
	};

#endif