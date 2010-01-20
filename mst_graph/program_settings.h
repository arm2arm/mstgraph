#ifndef PROGRAM_SETTINGS
#define PROGRAM_SETTINGS
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <exception>
#include <iostream>
struct programm_settings
{
    std::string m_snap_file;               // file name to analyse
    std::string m_path;                    // output path;
	int m_test_type; // is small test enabled?
	std::string m_HOP_file;

        struct TEstData{
            int dim;
            int ngb;
            int klink;
            TEstData(int d):dim(d){};
            std::string rhofile;
            std::string hsmlfile;
            std::string ngbfile;
            int flag;
            int smooth_flag;
        };
        int m_AnnTune;
//////////////////////////////////	
	struct TPathMask{
		std::string path; 
		std::string mask;
		std::string name;
		std::string get_Path(){return path;};
		std::string get_Name(){return path+mask+"_";};
		TPathMask(std::string n,std::string p, std::string m):name(n),path(p), mask(m){};
		TPathMask(){};
		};
	std::vector<std::string> m_snapvec;
	std::vector<TPathMask> m_modelvec;
        std::vector<TEstData> m_estvec;

	TPathMask m_HOPname, m_TESTname;
	////////////////////////////////
	programm_settings():m_test_type(false){};
    void load(const std::string &filename);
	void save(const std::string &filename){};
	int get_ModelCount(){return m_modelvec.size();};
	//std::string get_ModelName(){return m_modelvec[0].name;};
	void print(){
		BOOST_FOREACH(TPathMask &v,
			m_modelvec)
			std::cout<<
			"MODEL: "<<v.name<<
			"\nPATH: "<<v.path<<
			"\nMASK: "<<v.mask<<
			std::endl;
		};
	// Member functions for returning filenames: 
	std::string get_SNAPfile(int isnap,int imodel=0){
		return m_modelvec[imodel].get_Name()+m_snapvec[isnap];
		}
	std::string get_HOPfile(int isnap,int imodel=0)
		{
		return	m_HOPname.path+"/"+m_modelvec[imodel].mask+"_"+m_snapvec[isnap]+"."+m_HOPname.mask;
		}
};

extern  struct programm_settings pset;
#endif
