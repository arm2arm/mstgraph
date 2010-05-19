#ifndef _MSTGROUP_
#define _MSTGROUP_

#include <iostream>
#include <vector>
#include <cstring>
#include <stdlib.h>
using std::cout;
using std::endl;

////////////////////CGalaxy class ////////
class CMSTGroup{
public:

	CMSTGroup(){
		ResetAll();		
		//ID++;
		}
	void ResetAll()
		{		
		memset(com,0,sizeof(com));
		memset(wcom,0,sizeof(wcom));
		memset(vel,0,sizeof(vel));
		memset(wvel,0,sizeof(wvel));
		memset(npart,0,sizeof(npart));
		memset(Rxyz,0,sizeof(Rxyz));
		Ntotal=0;R90=0;R50=0;m_Rxx=0;
		wtotal=0.0;
		}
	/*void insert(int i, float x,float y,float z)
		{
		com[0]+=x;com[1]+=y;com[2]+=z;
		id.push_back(i);
		Ntotal=id.size();
		}*/
	void insert(int i, float x,float y,float z, float w)
		{		 
		com[0]+=x;com[1]+=y;com[2]+=z;
		wcom[0]+=x*w;
		wcom[1]+=y*w;
		wcom[2]+=z*w;
		wtotal+=w;
		//exit(0);
		id.push_back(i);
		Ntotal=id.size();
		}
	~CMSTGroup(){id.clear();};
	/// For sorting
	friend bool operator<(const CMSTGroup& left,const CMSTGroup& right);
	/*CMSTGroup& CMSTGroup::operator=(const CMSTGroup& p)
		{
		
		}*/
	friend std::ostream& operator<<(std::ostream& os, const CMSTGroup& g)
		{
		os.width(1);
		os<<g.ID;
		os.width(8);
		os<<" "<<g.Ntotal;

		os.precision(5);
		os.width(10);
		os.setf( std::ios::fixed, std::ios::floatfield ) ;
		os<<" "<<g.com[0];

		os.width(10);
		os<<" "<<g.com[1];

		os.width(10);
		os<<" "<<g.com[2];

		os.width(10);			
		os<<" "<<g.wcom[0];
		os<<" "<<g.wcom[1];
		os<<" "<<g.wcom[2];

		os.width(10);
		os<<" "<<g.wvel[0];
		os<<" "<<g.wvel[1];
		os<<" "<<g.wvel[2];


		os<<std::endl;

		return os;
		}
	friend std::istream &operator>>(std::istream& is,  CMSTGroup& g){
		
		is>>g.ID
		>>g.Ntotal
		>>g.com[0]
		>>g.com[1]
		>>g.com[2]

		>>g.wcom[0]
		>>g.wcom[1]
		>>g.wcom[2]

		>>g.wvel[0]
		>>g.wvel[1]
		>>g.wvel[2];


		return is;
	 
	 }
	void DoneInsert(int i=0)
		{
		ID=i;
		Ntotal=id.size();
		for(int j=0;j<3;j++)
		 {
		 com[j]/=double(Ntotal);
		 wcom[j]/=double(wtotal);
		 //cout<<wcom[j]<<" ";
		 //vel[j]/=double(Ntotal);
		 //wvel[j]/=(wtotal*Ntotal);
		 }
		//cout<<endl;
		}
	inline size_t size(){return id.size();};
	/////////////////////////////	
	std::vector<int> id;	
	double com[3], wcom[3], vel[3],wvel[3];
	float R90, R50;
	float m_Rxx;
	int ID;
	double wtotal;
	int npart[3];
	int Rxyz[3];
	int Ntotal;

	};

typedef std::vector<CMSTGroup> TMSTCat;

template <class T>
struct IfGt{
	T val;
	IfGt(const T n):val(n){};
	bool operator ()(const CMSTGroup& g){return g.Ntotal > val;};
	};
template <class T>
struct IfLt{
	T val;
	IfLt(const T n):val(n){};
	bool operator ()(const CMSTGroup& g){return g.Ntotal < val;};
	};
#endif


