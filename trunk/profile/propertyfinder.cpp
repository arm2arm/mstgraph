// profile.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//


#include "loader.h"
#include "options.h"
#include "utils.h"
#include "Logger.h"
#include "MSTree.h"
#include "PCA.h"
#include <numeric>
#include <valarray>
#include <functional>
#include "Sigma.h"
using std::string;
///////////////////////

int main(int argc, char* argv[])
	{
    scoped_timer timemme("Main program PropertyFinder :.....");
	COptions opt(argc, argv);

	if(opt.is_bad())
		return EXIT_FAILURE;
	unsigned int ParticleType=opt.m_type;
	typedef std::vector<CLoader*> TLoader;
	TLoader Lvec;
	cout<<"# snap Ngas<1 Ngas<2 MBHinit MBH  sigVstar sigVdisk"<<endl;
	for(size_t i=0;i<opt.m_snapshotList.size();i++)
		{
		unsigned int isnap=GetISnap(opt.m_snapshotList[i]);  
		if(!is_file_exist(opt.m_snapshotList[i])){cout<<"skipping...isnap="<<isnap<<endl;continue;}


		CLoader *pL=new CLoader(opt.m_snapshotList[i],6, true);//get vel as well
		CLoader *pLBH=new CLoader(opt.m_snapshotList[i],5, true);
		string path;
		GetSnapPath(opt.m_snapshotList[i], path);
		path=path+"/snap_m8_gal_sfr_bh_104";
		CLoader *pLBHINIT=new CLoader(path,5);
		
		
		///////////////////////////////
		//Do something on Snaps.
		float com[]={
		  pLBH->pPOS[0], pLBH->pPOS[1],pLBH->pPOS[2],
		  pLBH->pVEL[0],pLBH->pVEL[1],pLBH->pVEL[2]
		};
		pL->MoveToCOM(&com[0], 6);
		///////////////////////////////
		size_t ig1=0, ig2=0, istar=0;
		vector<int> itype(6,0);
		vector<double> sigV(6,0.0);
		for(size_t i=0;i<pL->size();i++)
		  {
		    
		    if(pL->pType[i]==0)
		      {
			if(pL->R(i)<1.0 )
			  {
			    ig1++;
			    
			  }
			if(pL->R(i)<2.0)
			  {
			    ig2++;
			  }
		      }
		    if(pL->pType[i]==4 && pL->R(i)<1.5)
		      {
			sigV[4]+=pL->Get3DVel(i);
			istar++;
			itype[4]++;
		      }
		    if(pL->pType[i]==2 && pL->R(i)<1.5)
		      {
			sigV[2]+=pL->Get3DVel(i);
			istar++;
			itype[2]++;
		      }
		  }
		for(size_t isig=0;isig<sigV.size();isig++)
		  if(itype[isig]>0)
		    sigV[isig]/=double(itype[isig]);
		///////////////////////////////
		cout.precision(10);
		int isnap=GetISnap(pL->m_fname);
		cout<<std::setw(12)<<std::fixed<<isnap<<" "
		    <<ig1<<" "<<ig2<<" "
		    <<pLBHINIT->GetBHMass()<<" "<<pLBH->GetBHMass()<<" "
		    <<sigV[2]<<" "<<sigV[4]<<endl;		
		///////////////////////////////
			{
			CSigma<double> sigma( &pL->pType[0],&pL->pPOS[0], &pL->pVEL[0],pL->size());
			sigma.m_fname="sigma_"+boost::lexical_cast<std::string>(isnap)+string("_4.txt");
			}
///////////////////////////////////

		delete pL;
		delete pLBH;
		delete pLBHINIT;
		}


	return 0;
	}

