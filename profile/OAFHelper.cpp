#include "OAFHelper.h"
#include <set>
using std::set;
/*COAFHelper::COAFHelper(void)
{
}
*/
COAFHelper::~COAFHelper(void)
	{
	///// Lets us dump the data.
	size_t Ns=snapu.size();
	Open('w');
	WriteHeader(idu_sel.size(),string(strOAF_52), 0, 0, 0.0, Ns);
	char *pData;
	//////////////////////////////////
	unsigned int icount=0,nbw=0;
	//////////////////////////////////
	unsigned int timeSize=Ns*sizeof(double);
	pData= new char[timeSize];
	vector<double> Time;
	for(id_unique_set::iterator itt=snapu.begin();itt!=snapu.end();++itt)
		{
		double number=(double)(*itt);
		memcpy((void*)&pData[icount], (void*)&(number), 8);
		Time.push_back((*itt));
		icount+=8;
		}

	WriteTime(pData, timeSize);

	delete [] pData;
	/////////////
	unsigned int DataSize=WriteDataStart(Ns), Nsc, iskip=0;
	icount=0;
	float var[NUM_RECS];
	pData= new char[4+52*Ns];
	memset(pData, 0, 4+52*Ns);
	set<T_indexID, setSnapLess> pset;
	for(id_unique_set::iterator it=idu_sel.begin();it!=idu_sel.end();++it)
		{
		int myID = (*it); 
		T_indexID ic0, ic1;
		boost::tuples::tie(ic0,ic1)=get<ID>(ptracks).equal_range(myID);

		if(ic0 ==ic1)
			{
			cout<<"Strange Error in ID sort"<<endl;
			exit(0);
			}

		Nsc=0;
		memmove((void*)&pData[0], (void*)&myID, 4);
		int count=0;
		/* MULTI THREADED THINGS */
		/*We need to sort them by time*/
		pset.clear();		  
		while( ic0 !=ic1 && ((*ic0).ID==myID))
			{
			pset.insert(ic0);
			ic0++;
			}

		set<T_indexID, setSnapLess>::iterator itSet=pset.begin();
		while( itSet !=pset.end()){
			////////////////////////////////
			//fPOS3:fVEL3:fMASS:fHSML:fRHO:fU:fZ2:fSFR
			var[0]=(*itSet)->data[X];
			var[1]=(*itSet)->data[Y];
			var[2]=(*itSet)->data[Z];

			var[3]=(*itSet)->data[VX];
			var[4]=(*itSet)->data[VY];
			var[5]=(*itSet)->data[VZ];

			var[6]=(*itSet)->data[MASS];
			var[7]=(*itSet)->data[HSML];
			var[8]=(*itSet)->data[RHO];
			var[9]=(*itSet)->data[U];
			var[10]=(*itSet)->data[Zmg];
			var[11]=(*itSet)->data[Zms];
			var[12]=(*itSet)->data[SFR];
			var[13]=(*itSet)->data[MASS];


			if(Nsc > Ns*52 )
				{
				cerr<<"ERROR: Nsc > Ns*52 : "<<Nsc<<endl;
				cerr<<"ID: "<<myID<<" Nsc= "<<Nsc<<" we got:"<<Nsc/32<<" should be  Ns= "<<Ns
					<<" SNAP: "<<(*itSet)->snap<<endl;
				cerr<<(*(*itSet));
				//error_flag=1;
				iskip++;
				break;
				}else
				{

				memmove((void*)&pData[4+Nsc], 
					(void*)&var[0], 52);	  
				Nsc+=32;
					}
				////////////////////////////////
				++itSet;
			}
		WriteOneParticle(pData, 4+52*Ns);
		nbw+=(4+52*Ns);
		///////////////////////////////////////
		if(icount%1000 ==0 )//&&VERBOSE_flag)
			{
			cout<<icount<<" ) ";
			cout<<"OAF Dumping the orbit for ID: "<< myID<<endl;
			}

		///////////////////////////////////////
		icount++;

		}

	WriteDataStop();
	Close();
	delete [] pData;
	std::cout.precision(10);
	cout<<"####################"<<endl;
	cout<<"Dumped n= "<<icount<<" orbits"<<endl;
	cout<<"DATA block total size :"<<nbw/(1024.0*1024.0)<<"Mb"<<endl;
	cout<<"Number of  BAD particles with duplicate IDS: "<<iskip<<endl;
	cout<<"####################"<<endl;
	}
bool COAFHelper::is_id_inselection(int id)
	{
	if(idu_sel.find(id)==idu_sel.end())
		return false;
	return true;
	}
bool COAFHelper::insert(int id, int snap, unsigned char type, float* pos, float *vel,
						float mass, float hsml,float rho, float u, float *z, float sfr)
	{
	bool is_id_there=is_id_inselection(id);
	if(is_id_there)
		{
		snapu.insert(snap);
		ptracks.insert(
			strParticle(id, snap, type, pos, vel,
			mass,hsml,rho,u,z,sfr)
			); 
		}
	return is_id_there;
	}



