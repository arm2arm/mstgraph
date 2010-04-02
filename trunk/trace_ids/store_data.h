#ifndef _MY_STOREDATA_
#define _MY_STOREDATA_


#include <fstream> 
#include <iostream> 
#include <boost/iostreams/copy.hpp> 
#include <boost/iostreams/filter/bzip2.hpp> 
#include <boost/iostreams/device/file.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/combine.hpp> 
#include "loader.h"

class CStoreData{
 private:
  std::string  m_fname;
 public:
	~CStoreData(){std::cout<<"...DONE"<<std::endl;};
	CStoreData(){};
	CStoreData(std::string filename,CLoader* L1, CLoader* L2, particlesID_set &ids_inboth){
			DumpIntersection(filename,L1,L2, ids_inboth);
			};
	void DumpIntersection(string filename,CLoader* L1, CLoader* L2, particlesID_set &ids_inboth)
		{
		m_fname=filename+".bz2";
		std::cout<<"\n Writing: "<<m_fname;
		namespace BI = boost::iostreams;

		std::fstream myFile(m_fname.c_str(), std::ios::binary|std::ios::out);
		BI::filtering_stream<BI::output> my_filter; 
		my_filter.push(BI::bzip2_compressor()) ; 
		my_filter.push(myFile) ; 

		my_filter<<L1->m_fname<<"\t"<<L1->m_nelem<<"\t"<<L1->m_ptype<<endl;
		my_filter<<L2->m_fname<<"\t"<<L2->m_nelem<<"\t"<<L2->m_ptype<<endl;

		boost::multi_index::index<particlesID_set,ID>::type& L1_ID_index=
			get<ID>(L1->data);
		boost::multi_index::index<particlesID_set,ID>::type& L2_ID_index=
			get<ID>(L2->data);
		boost::multi_index::index<particlesID_set,ID>::type::iterator  L1_ID_index_it;
		boost::multi_index::index<particlesID_set,ID>::type::iterator  L2_ID_index_it;
		particlesID_set::iterator it=ids_inboth.begin();
		while(it!=ids_inboth.end())
			{
			L1_ID_index_it=L1_ID_index.find(it->ID);
			L2_ID_index_it=L2_ID_index.find(it->ID);
			my_filter<<it->ID<<" "<<L1_ID_index_it->IDf<<" "<<L2_ID_index_it->IDf<<std::endl;
			it++;
			}
		}
	void test_bzip()
		{
		namespace BI = boost::iostreams;
			{
			string fname="test.bz2";
				{
				  std::fstream myfile(fname.c_str(), std::ios::binary|std::ios::out); 
				BI::filtering_stream<BI::output> my_filter; 
				my_filter.push(BI::bzip2_compressor()) ; 
				//				my_filter.push(std::fstream(fname.c_str(), std::ios::binary|std::ios::out)) ; 
				my_filter.push(myfile);
				my_filter << "test";
				}
			}
		}
	};

#endif
