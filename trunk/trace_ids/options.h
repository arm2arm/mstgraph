#ifndef _MYOPTIONS_
#define _MYOPTIONS_
#include <iostream>
#include <exception>
#include <fstream>
#include <cstdlib>
#include <map>
#include <set>
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ostream_iterator;
using std::map;
using std::set;

class COptions{
public:
	COptions(int argc, char* argv[]);
	const int is_bad(){return m_status;};
//////////////////////////////////////////
	vector<string> m_snapshotList;
	set<int> m_IDlist;
	string m_file_out;
	string m_file_IDlist;
	int m_type;
private:
	void ReadID(std::string file);
	int m_status;
	};

#endif