//=======================================================================
// Copyright 2009 Arman Khalatyan, OAMP/LAM  
// distributed under GPL license
//=======================================================================
#include "mst_graph.h"
#include "kdtree2.hpp"
#include "functions.h"
#include "Render.h"
#include "HOP.h"
#include "GetEst.h"
/////////////////////////////
//#define ND_GROUPS 2
#define LOAD_DEBUG 0
#define SMALL_TEST 0
#define MIN_NGRP 50
#define MAXNGB 64.0 
std::string base;//this is used for outputs
#define RAND_FRACTION 1.2
///////////////////////////////
//
// define, for convenience a 2d array of floats. 
//  
typedef boost::multi_array<MyFloat,2> array2dfloat;

///////////////////////////////
#include "readers.h"
void read_gadget(string fname)
	{
	All.ONLY_TYPE=4;
	read_ic12(fname.c_str(), true);
/*	MyFloat COM[]={0.0,0.0,0.0}, totEst=0.0;
	for(int i=0;i<All.NumPart;i++)
		{
		COM[0]+=Part[i].Pos[0]*Part[i].Est;		
		COM[1]+=Part[i].Pos[1]*Part[i].Est;		
		COM[2]+=Part[i].Pos[2]*Part[i].Est;		
		totEst+=Part[i].Est;
		}
	COM[0]/=totEst;
	COM[1]/=totEst;
	COM[2]/=totEst;
	for(int i=0;i<All.NumPart;i++)
		{
		Part[i].Pos[0]-=(float)COM[0];
		Part[i].Pos[1]-=(float)COM[1];
		Part[i].Pos[2]-=(float)COM[2];
		}
		*/
	}

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <exception>
struct programm_settings
{
    std::string m_snap_file;               // file name to analyse
    std::string m_path;                    // output path;
	int m_test_type; // is small test enabled?
	std::string m_HOP_file;

//////////////////////////////////
	struct TPathMask{
		std::string path; 
		std::string mask;
		std::string name;
		TPathMask(std::string n,std::string p, std::string m):name(n),path(p), mask(m){};
		};
	std::set<std::string> m_snapvec;
	std::vector<TPathMask> m_modelvec;
	////////////////////////////////
	programm_settings():m_test_type(false){};
    void load(const std::string &filename);
	void save(const std::string &filename){};
	void print(){
		BOOST_FOREACH(TPathMask &v,
			m_modelvec)
			cout<<
			"MODEL: "<<v.name<<
			"\nPATH: "<<v.path<<
			"\nMASK: "<<v.mask<<
			endl;
		};
} pset;
void programm_settings::load(const std::string &filename)
	{
	using boost::property_tree::ptree;
	ptree pt;
	std::string ext(filename.end()-3,filename.end());
	if(ext=="ini")
		{
		read_ini(filename, pt);
		m_test_type = pt.get<int>("OPTIONS.test_flag");
		if(m_test_type==1)
			m_snap_file = pt.get<std::string>("OPTIONS.test_file");
		else
			m_snap_file = pt.get<std::string>("OPTIONS.file");

		m_path = pt.get<std::string>("OPTIONS.path");
		m_HOP_file = pt.get<std::string>("OPTIONS.HOP_file");
		}
	else
		{
		
		read_xml(filename, pt);
		int test_flag = pt.get("OPTIONS.test_flag", 0);

		BOOST_FOREACH(ptree::value_type &v,
			pt.get_child("OPTIONS.snaps"))
			m_snapvec.insert(v.second.data());
	
		BOOST_FOREACH(ptree::value_type &m,
			pt.get_child("OPTIONS.models"))
			{
			m_modelvec.push_back(TPathMask(m.second.data(),"",""));
			}
		if(m_modelvec.size() > 0)
			{
			int imodel=0;
			BOOST_FOREACH(ptree::value_type &p,
				pt.get_child("OPTIONS.paths"))
				{
				m_modelvec[imodel++].path=p.second.data();
				}
			imodel=0;
			BOOST_FOREACH(ptree::value_type &p,
				pt.get_child("OPTIONS.masks"))
				{
				m_modelvec[imodel++].mask=p.second.data();
				}

			}
		
		}
	print();
	}

void myreadini(string fini)
	{
	try
		{
		pset.load(fini);
		//pset.save(fini);
		std::cout << "Success\n";
		}
	catch (std::exception &e)
		{
		std::cout << "Error: " << e.what() << "\n";
		}
	};
void myreadxml(string fxml)
	{
	try
		{
		pset.load(fxml);
		//pset.save(fini);
		std::cout << "Success\n";
		}
	catch (std::exception &e)
		{
		std::cout << "Error: " << e.what() << "\n";
		}
	};

///////////////////////////////////////////////////////////////////
//"snap_gal_sfr_0450.ascii"	
int main(int argc,char **argv) {

	if(argc==2)
		myreadini(argv[1]);
	// Read and check command line parameters.
//	cimg_usage("Compute a HOP over the particles with given Est and Rho files");
//	const char *file_ascii = cimg_option("-ascii_file",(char*)0,"Input in Ascii format");
//	const char *file_gadget2;
//	if( SMALL_TEST)
//		file_gadget2= cimg_option("-gadget2_file","C:\\arm2arm\\DATA\\MODEL7\\MODELS\\MODEL7\\RUNG2\\SNAPS\\test_0450","Input in Gadget-2 format");
//	else
//		file_gadget2= cimg_option("-gadget2_file","C:\\arm2arm\\DATA\\MODEL7\\MODELS\\MODEL7\\RUNG2\\SNAPS\\snap_gal_sfr_0450","Input in Gadget-2 format");
//	const char *file_base  = cimg_option("-base4files","C:/arm2arm/DATA/MODEL7/MST_GRAPH/","Base for input and output");
//	const char *HOP_file  = cimg_option("-hopfile","hop_0450","HOP algorithm output,eg: groups and their IDs.");
/*	const bool visu     = cimg_option("-visu",true,"Visualization mode");
	const int shape  = cimg_option("-s",1,"shape [0,6]");
	const int profile  = cimg_option("-p",0,"profile [0,7]");
*/
	////////////////////////////////////////////////////////////////////

	base=pset.m_path;
	string hopfile=base+string("/")+pset.m_HOP_file;
		
	read_gadget(pset.m_snap_file);

	////////////////////////////////////////////////
	/// Smooth Est with 64 Ngb
	GetEst est_me;
	est_me.Run_SPHEst();
	exit(0);
	//SmoothSph();
	/////////////////////// GET setAB This is the HOP stuff based on Enbid Density  one can test also for RHO by SPH//////////////////////////////////////
	//exit(0);
	MyFloat alpha=0.05;
	CHOP *hop=new CHOP(alpha, MIN_NGRP);

	int seeds=hop->get_seeds();
	cout<<"We got "<<seeds<<" seeds out of "<< All.NumPart<<endl;


	hop->write_catalogues(hopfile);

	///////////////////////////////////////////////////////////////
	// Here we need to cut and store the catalogues
	///////////////////////////////////////////////////////


	////////////////////////////////
	return EXIT_SUCCESS;
	}
