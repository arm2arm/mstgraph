#ifndef _MSTGROUP_
#define _MSTGROUP_

#include <iostream>
#include <vector>
#include <cstring>


////////////////////CGalaxy class ////////
class CMSTGroup{
public:
	CMSTGroup(){
		ResetAll();		
		};
	CMSTGroup(int comp, typeVecData *datain):data(datain){
		ResetAll();
		ID=comp;
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
	void CompileMyProperties(float *pPOS, float *pVEL, float *pW, int np)
		{
		for(int i=0,pc=0;i<Ntotal;i++)
			{
			pc=id[i]*3;
			double w=pW[id[i]]*pW[id[i]];
			for(int j=0;j<3;j++)
				{
				com[j]+=pPOS[pc+j];
				wcom[j]+=(pPOS[pc+j]*w);
				vel[j]+=(pVEL[pc+j]);
				wvel[j]+=(pVEL[pc+j]*w);
				}
			wtotal+=w;
			}
		for(int j=0;j<3;j++)
			{
			com[j]/=double(Ntotal);
			wcom[j]/=(double(Ntotal)*wtotal);
			vel[j]/=double(Ntotal);
			wvel[j]/=(double(Ntotal)*wtotal);
			}
		};
	~CMSTGroup(){};

 friend std::ostream& operator<<(std::ostream& os, const CMSTGroup& g)
    {
      os.width(1);
  	  os<<g.ID;
      os.width(8);
      os<<" "<<g.Ntotal;
      
	  os.precision(4);
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
 void Insert(int i){
	 id.push_back(i);
	
	 double w=(*data)[i].w;
	 w*=w;
	 for(int j=0;j<3;j++)
		 {
		 com[j]+=(*data)[i].pos[j];
		 wcom[j]+=((*data)[i].pos[j]*w);
		 vel[j]+=((*data)[i].vel[j]);
		 wvel[j]+=((*data)[i].vel[j]*w);
		 }
	 wtotal+=w;

	 }
 void DoneInsert()
	 {
	 Ntotal=id.size();
	 for(int j=0;j<3;j++)
		 {
		 com[j]/=double(Ntotal);
		 wcom[j]/=(wtotal*Ntotal);
		 vel[j]/=double(Ntotal);
		 wvel[j]/=(wtotal*Ntotal);
		 }
	
	 }
	/////////////////////////////	
	vector<int> id;
	typeVecData *data;
	double com[3], wcom[3], vel[3],wvel[3];
	float R90, R50;
	float m_Rxx;
	int ID;
	double wtotal;
	int npart[3];
	int Rxyz[3];
	int Ntotal;
	};


#endif


