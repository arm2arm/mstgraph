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
	std::vector<string> m_snapshotList;
	std::vector<std::string>  m_ApFilelist;
	std::vector<std::vector<int> >  m_IDlistvec;
	std::vector<int>  m_IDlist;
	string m_file_out;
	std::vector<std::string> m_IDfilelist;
	int m_type;
	float m_eps;
	int  m_min_npart;
	int m_NGB;
	float m_Rmax;
	bool m_updatelog;
	bool m_OAF;
private:
	bool ParseSnapshotLists(std::vector<string> &strinout);
	void ReadID(std::vector<std::string > &filevec);
	int m_status;

	typedef std::vector<string> Tstrvec;
	typedef Tstrvec::iterator TstrvecIT;
	};

#endif


