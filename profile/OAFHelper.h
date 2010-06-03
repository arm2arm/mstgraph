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

enum TAGS{X, Y, Z, VX, VY, VZ, TYPE,MASS, HSML, RHO, U, Zmg,Zms, SFR, NUM_RECS};//Zm =Z in  gadget, gas and stars
struct strParticle{
	strParticle(int id, int snap,unsigned char type, float* pos, float *vel,
		float mass, float hsml,float rho, float u, float *z, float sfr):ID(id), 
		snap(snap), type(type)
		{
		data.resize(NUM_RECS);
		data[X]=pos[0];
		data[Y]=pos[1];
		data[Z]=pos[2];
		data[VX]=vel[0];
		data[VY]=vel[1];
		data[VZ]=vel[2];
		data[MASS]=mass;
		data[TYPE]=type;
		data[HSML]=hsml;
		data[RHO]=rho;
		data[U]=u;
		data[Zmg]=z[0];
		data[Zms]=z[1];
		data[SFR]=sfr;
		};// 3Pos, 3Vel, MASS, HSML, RHO, U, Z, SFR

	int    ID;
	unsigned short snap;
	std::vector<float> data;
inline 	float  GetR(void) const{return std::sqrt(data[0]*data[0]+data[1]*data[1]+data[2]*data[2]);};
inline 	float  R(void) const{return GetR();};

	size_t size(){return data.size();};
	char type;

    friend std::ostream& operator<<(std::ostream& os, const strParticle& e)
    {
      os<<"SNAP: "<<e.snap<<" ID: "<<e.ID<<" R: "<<e.R()<<std::endl;
     /* os<<e.snap;
      os.width(16);
      os<<"\t"<<e.id;
      os.precision(4);
      os.width(16);
      os.setf( std::ios::fixed, std::ios::floatfield ) ;
      os<<" "<<e.Pos[0];
      os.width(16);
      os<<"\t"<<e.Pos[1];
      os.width(16);
      os<<"\t"<<e.Pos[2];
      os.width(10);			
      os<<"\t"<<e.GetR()<<std::endl;
*/
      return os;
    }
 
	};

class COAFHelper:public COAFFile
	{
	struct ID{};struct SNAP{};struct R{};
	public:
		COAFHelper(std::string file="oaf"):COAFFile(file){};
		~COAFHelper();
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
		id_unique_set idu_sel, snapu;
		particle_set ptracks;

/////////////// USED TO SORT THE MULTITHREADED SETS ////
		typedef boost::multi_index::index<particle_set,ID>::type::iterator T_indexID;
		struct setSnapLess {
			bool operator( )(const T_indexID p1,
				const T_indexID p2) {
					return( (*p1).snap < (*p2).snap);
				}
			};

	};

#endif


