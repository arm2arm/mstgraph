
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

using namespace boost::program_options;
using boost::lexical_cast;
using boost::bad_lexical_cast;

#include "options.h"

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


COptions::COptions(int argc, char* argv[]):m_status(0)
	{
	m_status=EXIT_SUCCESS;
	try {
		namespace po = boost::program_options;
		//typedef string TIDlist;
		typedef int TIDlist;
		std::vector<TIDlist> IDlist;
		po::options_description commands("Allowed options");
		commands.add_options()
			("snapshotList", po::value< vector<string> >(&m_snapshotList)->multitoken(), "Which snapshots to trace: ex. --snapshotList=snap_001 snap_002 snap_003 ")
			("IDlist",       po::value< vector<TIDlist> >(&IDlist)->multitoken(), "Which IDs to trace: ex. --IDlist=0 1 200 -340")
			("out-file", po::value< string>(&m_file_out)->default_value(string("intersect.idx")), "out file")
			("type", po::value< int>(&m_type)->default_value(4), "particle type to track")
			("help","print help")
			;

		variables_map vm;
		store(parse_command_line(argc, argv, commands, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
		//		store(command_line_parser(argc, argv).options(commands).run(), vm);// this is the simplyfied interface of the parser
		notify(vm);

	//	printInputs(vm);

		if (vm.count("help")) {
			cout << commands << "\n";
			m_status=EXIT_FAILURE;
			}

		if(vm.count("snapshotList")) {
			//ostream_iterator<string> it(cout, "\n ");
			cout << "snapshot list detected" << endl;
			cout << "We will trace following snapshots:\n ";
			for(std::vector<string>::iterator it=m_snapshotList.begin();it<m_snapshotList.end(); it++)
				{
				cout<<(*it);
				if(is_file_exist(*it))
					cout<<"..exist."<<endl;
				else{
					cout<<"does not exist."<<endl;
					m_status=EXIT_FAILURE;
					}
				}
			cout << endl << endl;
			if(m_snapshotList.size()<2)
				{
				cout<<"Error: The minimum number of snapshots is 2.\nPlease enter at least 2 different snapshots."<<endl;
				m_status=EXIT_FAILURE;
				}
			}
		if(vm.count("IDlist")) {
			ostream_iterator<TIDlist> it(cout, ", ");
			cout << "ID list detected with:" << endl;
			//cout << "These IDs we will trace: ";
			BOOST_FOREACH(TIDlist val, IDlist)
				{
				m_IDlist.push_back(boost::lexical_cast<int>(val));
				};
			cout<<"numID="<<m_IDlist.size()<<" entries.";
			cout << endl << endl;
			}
		}
	catch (const std::exception& e) {
		cout << e.what() << "\n";
		m_status=EXIT_FAILURE;
		}

	}

/*void COptions::printInputs(variables_map vm)
{
int i=0;
variables_map::iterator it=vm.begin();
cout << "\n***List of commands detected: " << endl;
for(it;it!=vm.end();++it) {
cout <<"Command[" << i << "]: " << it->first << endl;
++i;
}
cout << endl;
}
*/