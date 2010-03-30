#ifndef _PARTICLE_ID_
#define _PARTICLE_ID_
#pragma warning(disable : 4996)
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <algorithm>
#include <iostream>
#include <ostream>
#include <iterator>
#include <string>

using boost::multi_index_container;
using namespace boost::multi_index;

struct particleID
	{
	int           ID;// real ID for particle from Gadget2 file "ID" block
	unsigned int  IDf;// postition in the file 
	particleID(int id,const unsigned int idf):ID(id),IDf(idf){}
	bool operator<(const particleID& p)const { return ID<p.ID;}
	unsigned int getByGasID()const {return (ID&0x0FFF);};
	friend std::ostream & operator<<(std::ostream& , const particleID& p );
	};

std::ostream & operator << (std::ostream &os, const particleID &p)
	{
	os <<"IDf: "<< p.IDf <<" < == > ID: "<< p.ID <<std::endl;
	return os;
	}

struct ID{};
struct IDf{};
struct IDgas{};

typedef multi_index_container<
	particleID,
	indexed_by<
		ordered_unique<
			tag<IDf>,  BOOST_MULTI_INDEX_MEMBER(particleID,unsigned int,IDf)>,
		ordered_non_unique<
			tag<ID>,BOOST_MULTI_INDEX_MEMBER(particleID,int,ID)>,
		ordered_non_unique<
			tag<IDgas>,BOOST_MULTI_INDEX_CONST_MEM_FUN(particleID,unsigned int,getByGasID)> 
	>
> particlesID_set;


template<typename T> struct strComparator;
template<> 
struct strComparator<ID>{
	bool operator () (const particleID& id1, const particleID& id2) const
		{
		return id1.ID<id2.ID;
		}
	};
template<> 
struct strComparator<IDf>{	
	bool operator () (const particleID& id1, const particleID& id2) const
		{
		return id1.IDf<id2.IDf; 
		}

	};
template<> 
struct strComparator<IDgas>{	
	bool operator () (const particleID& id1, const particleID& id2) const
		{
		return id1.getByGasID()<id2.getByGasID(); 
		}

	};

template<typename Tag,typename MultiIndexContainer>
void print_out_by(
				  const MultiIndexContainer& s,
				  Tag* =0 /* fixes a MSVC++ 6.0 bug with implicit template function parms */
				  )
	{
	/* obtain a reference to the index tagged by Tag */

	const typename boost::multi_index::index<MultiIndexContainer,Tag>::type& i=
		get<Tag>(s);

	typedef typename MultiIndexContainer::value_type value_type;

	/* dump the elements of the index to cout */

	std::copy(i.begin(),i.end(),std::ostream_iterator<value_type>(std::cout));
	}
#endif