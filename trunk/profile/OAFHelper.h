#ifndef _OAFHELPER_
#define _OAFHELPER_
#include "oaffile.h"
#include <vector>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>

namespace MI = boost::multi_index;
using boost::multi_index_container;
using namespace boost::multi_index;


struct strParticle{
	enum TAGS{X, Y, Z, VX, VY, VZ, MASS, HSML, RHO, U, Zm, SFR};//Zm =Z in  gadget
	strParticle(int id)
		{data.resize(3+3+3+3);};// 3Pos, 3Vel, MASS, HSML, RHO, U, Z, SFR
	int    ID;
	unsigned short snap;
	std::vector<float> data;
	float  GetR(void) const{return std::sqrt(data[0]*data[0]+data[1]*data[1]+data[2]*data[2]);};
	size_t size(){return data.size();};
	char type;
	};

class COAFHelper:public COAFFile
	{
	struct ID{};struct SNAP{};struct R{};
	public:
		COAFHelper(std::string file="oaf"):COAFFile(file){};
		~COAFHelper(){};
		bool insert(int id, int snap,unsigned char type, float* pos, float *vel,
			float mass, float hsml,float rho, float u, float *z, float sfr);
		void SetID(std::vector<int> &id)
			{
			std::copy(
				id.begin(),
				id.end(),
				std::inserter(idu_sel, idu_sel.end())
				);

			};
		bool is_id_inselection(int id);
	private:
		typedef MI::multi_index_container<
			strParticle,
			MI::indexed_by<
			MI::ordered_non_unique<
			tag<ID>,  BOOST_MULTI_INDEX_MEMBER(strParticle,int,ID)>,
			MI::ordered_non_unique<
			tag<SNAP>,BOOST_MULTI_INDEX_MEMBER(strParticle,unsigned short,snap)>
#ifdef R_INDEX
			,MI::ordered_non_unique<

			tag<R>, BOOST_MULTI_INDEX_MEMBER(strParticle,float,R)>
#endif
			>

		> particle_set;
		////////////////////////////////////////////////
		typedef MI::multi_index_container<
			int,
			MI::indexed_by<
			MI::ordered_unique<MI::identity<int> >
			>
		> id_unique_set;
		id_unique_set idu_sel;
	};

#endif