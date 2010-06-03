
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
#include "utils.h"

#include <fstream>
#include <string>
#include <set>

#include <boost/algorithm/string.hpp>// for split
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>  // also functional,iostream,iterator,string
namespace bfs = boost::filesystem;


void COptions::ReadID(std::vector< std::string > &filevec)
	{
	
	}

COptions::COptions(int argc, char* argv[]):m_status(0)
	{
	m_status=EXIT_FAILURE;
	try {
		namespace po = boost::program_options;
		//typedef string TIDlist;
		typedef int TIDlist;
		std::vector<TIDlist> IDlist;
		po::options_description commands("Allowed options");
		commands.add_options()
			("snapshotList", po::value< vector<string> >(&m_snapshotList)->multitoken(), "Which snapshots to profile: ex. --snapshotList=snap_001 snap_002 snap_003 ")
			("IDlist",       po::value< vector<TIDlist> >(&IDlist)->multitoken(), "Which IDs to trace: ex. --IDlist=0 1 200 -340")
			("IDfilelist",      po::value< vector<string> >(&m_IDfilelist)->multitoken(), "Which IDs to trace from file: ex. --IDfilelist=idFile.txt, where file has a integers inside the file")
			("out-file", po::value< string>(&m_file_out)->default_value(string("dump.log")), "out file, the result file if it exist the code will update if update=true(default value)")
			("update", po::value< bool>(&m_updatelog)->default_value(1), " update out file 0/1= No/Yes, if =0 then it will overwrite  out file")
			("type", po::value< int>(&m_type)->default_value(6), "particle type to track, if 6 then all particles")
			("r-max,r", po::value< float>(&m_Rmax)->default_value(20.0f, "20"), "MSTFOF linking length ")
			("ngb,n", po::value< int>(&m_NGB)->default_value(15), "NGB number for density estimation ")
			("eps,e", po::value< float>(&m_eps)->default_value(0.05f, "0.05"), "MSTFOF linking length ")
			("min-npart,m", po::value< int>(&m_min_npart)->default_value(100), "MSTFOF group minimum length ")			
			("mst-file,o", po::value< string>(&m_file_out)->default_value(string("mstfof")), "MSTFOF group out file")
			("OAF", po::value< bool>(&m_OAF)->default_value(0), " dump OAF trace file for the particles")
			("help,h","print help")
			;

		variables_map vm;
		store(parse_command_line(argc, argv, commands, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
		//		store(command_line_parser(argc, argv).options(commands).run(), vm);// this is the simplyfied interface of the parser
		notify(vm);

	//	printInputs(vm);

		if(vm.count("snapshotList")) {
			cout << "# snapshot list detected" << endl;
			cout<<"# Parsing snapshot lists..."<<endl;
			ParseSnapshotLists(m_snapshotList);
			cout << "# We will trace following snapshots:\n";
			for(std::vector<string>::iterator it=m_snapshotList.begin();it<m_snapshotList.end(); it++)
				{
				  cout<<"# "<<(*it);
				if(is_file_exist(*it)){
					cout<<"..exist."<<endl;
					m_status=EXIT_SUCCESS;
					}
				else{
					cout<<"does NOT exist."<<endl;
					//m_status=EXIT_FAILURE;
					}
				}
			//cout << endl << endl;
			if(m_snapshotList.size()<1)
				{
				cout<<"Error: The minimum number of snapshots is 1."<<endl;
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
		if(vm.count("IDfilelist")) {
			//ParseSnapshotLists(m_IDfilelist);
			std::set<std::string> result;
			for(TstrvecIT it=m_IDfilelist.begin();it<m_IDfilelist.end(); it++)
				boost::split(  result,*it,boost::is_any_of("\t ") );
			//count the entries
			m_IDfilelist.clear();
			std::unique_copy(result.begin(),result.end(),std::back_inserter(m_IDfilelist));

			
			ReadID(m_IDfilelist);
			cout<<"numID="<<m_IDlistvec.size()<<" entries.";
			cout << endl << endl;
			}
		if (vm.count("help")) {
			cout << std::fixed << std::setprecision( 3 )<<commands << "\n";
			m_status=EXIT_FAILURE;
			}
		if(m_status==EXIT_FAILURE)
			cout << std::fixed << std::setprecision( 3 )<<commands << "\n";
		}
	catch (const std::exception& e) {
		cout << e.what() << "\n";
		m_status=EXIT_FAILURE;
		}
	
	}



struct match : public std::unary_function<bfs::directory_entry,bool> {
	match(const std::string pat="snap_[0-2][0-9][0-9][0-9]"):pat(pat){}
	std::string pat;
    bool operator()(const bfs::directory_entry& d) const {
        std::string fn(d.filename());
        return boost::regex_match(fn.begin(), fn.end(), boost::regex(pat));
    }
};
template <class Tstrvec>
class CFileNameParser
{
public:
 CFileNameParser(std::string& path, Tstrvec& outlist)
  : m_path(path)
  , m_outlist(outlist)
 {}

 void operator()(const std::string& str) const
 {
  std::string fname = m_path + "/" + str; 
  boost::algorithm::replace_all(fname, "\\", "/"); 
  if( std::find(m_outlist.begin(),m_outlist.end(), fname) == m_outlist.end()) 
	  if(!fname.empty())m_outlist.insert(fname);  
 }

private:
 std::string m_path;
 Tstrvec& m_outlist;
};


bool COptions::ParseSnapshotLists(std::vector<string> &strinout){
	bool state=true;
	
	std::set<std::string> result;
	std::string pat, path;
	for(TstrvecIT it=strinout.begin();it<strinout.end(); it++)
		boost::split(  result,*it,boost::is_any_of("\t ") );
	//count the entries
	strinout.clear();
	std::unique_copy(result.begin(),result.end(),std::back_inserter(strinout));
	result.clear();
	for(TstrvecIT it=strinout.begin();it<strinout.end(); it++){
		if(is_file_exist(*it)){
			boost::algorithm::replace_all((*it), "\\", "/");
			if(std::find(result.begin(),result.end(),*it)==result.end())
				result.insert(*it);
			}else{
			//get path
			if((*it).empty())continue;
			bfs::path full_path( *it );
			path=full_path.branch_path().string();		
			pat=full_path.filename();
			boost::algorithm::replace_all(pat, "?", "[0-9]");
			std::list<std::string> partial_result;
			match PatternMatch(pat);
			transform_if(bfs::directory_iterator(path), bfs::directory_iterator(),
				std::back_inserter(partial_result),
				PatternMatch,
				mem_fun_ref(&bfs::directory_entry::filename));
			std::for_each(partial_result.begin(), 
				partial_result.end(), 
				CFileNameParser<std::set<std::string> >(path, result) );
			}
		};
	//std::copy(result.begin(), result.end(),std::ostream_iterator<std::string>(std::cout, "\n"));
	strinout.clear();
	std::unique_copy(result.begin(),result.end(),std::back_inserter(strinout));
	return state;
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
