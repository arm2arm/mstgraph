// trace_ids.cpp : Defines the entry point for the console application.
//

#include "loader.h"
#include "options.h"
#include "utils.h"
#include "store_data.h"

template<typename Tag,typename MultiIndexContainer>
void intersect_by(
				  const MultiIndexContainer& L1,const MultiIndexContainer& L2,MultiIndexContainer &s/*result*/, std::string msg="",
				  Tag* =0 /* fixes a MSVC++ 6.0 bug with implicit template function parms */
				  )
	{
	/* obtain a reference to the index tagged by Tag */

	const typename boost::multi_index::index<MultiIndexContainer,Tag>::type& L1_ID_index=
		get<Tag>(L1);
	const typename boost::multi_index::index<MultiIndexContainer,Tag>::type& L2_ID_index=
		get<Tag>(L2);
	typename boost::multi_index::index<MultiIndexContainer,Tag>::type::iterator  L1_ID_index_it;
	typename boost::multi_index::index<MultiIndexContainer,Tag>::type::iterator  L2_ID_index_it;
	//particlesID_set s;
	s.clear();
	std::set_intersection(
		L1_ID_index.begin(), 
		L1_ID_index.end(), 
		L2_ID_index.begin(), 
		L2_ID_index.end(),
		std::inserter(s, s.begin()), strComparator<Tag>()
		);

	
	}

int main(int argc, char* argv[])
	{
	COptions opt(argc, argv);
	if(opt.is_bad())
		return EXIT_FAILURE;
	unsigned int ParticleType=opt.m_type;
	typedef std::vector<CLoader*> TLoader;
	TLoader Lvec;

	for(size_t i=0;i<opt.m_snapshotList.size();i++)
		{
		CLoader *pL=new CLoader(opt.m_snapshotList[i], ParticleType);
		(*pL).PrintStat();
		Lvec.push_back(pL);
		}

		{
		for(size_t i=0;i<Lvec.size()-1;i++)
			{
			CLoader* L1, *L2;
			L1=Lvec[i];L2=Lvec[i+1];
			particlesID_set ids_inboth;
			intersect_by<ID>(L1->data, L2->data,ids_inboth);
			std::cout<<"by ID: We got " <<ids_inboth.size()<<
				" out of L1.size= "<<
				L1->m_nelem<<
				" and L2.size= "<<
				L2->m_nelem<<std::endl;

			/*intersect_by<IDgas>(L1->data, L2->data,ids_inboth);
			std::cout<<"by IDgas We got: " <<ids_inboth.size()<<
				" out of L1.size= "<<
				L1->m_nelem<<
				" and L2.size= "<<
				L2->m_nelem<<std::endl;
			*/
			CStoreData store_now(opt.m_file_out+string("_byIDstar"), L1, L2, ids_inboth);
			}
			
		
		}
		/// Lets free the memory
		delete_pointers_list< TLoader >(Lvec);
			
		//print_out_by<IDf>(L1.data);
		//intersect_by<IDf>(L1.data, L2.data, "By IDf: ");
		return 0;
	}

