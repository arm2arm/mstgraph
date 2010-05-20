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

	for(size_t i=0;i<opt.m_snapshotList.size();i++)
		{
		unsigned int isnap=GetISnap(opt.m_snapshotList[i]);  
		if(!is_file_exist(opt.m_snapshotList[i])){cout<<"skipping...isnap="<<isnap<<endl;continue;}


		CLoader *pL=new CLoader(opt.m_snapshotList[i],0);
		CLoader *pLBH=new CLoader(opt.m_snapshotList[i],5);

		///////////////////////////////
		//Do something on Snaps.
		pL->MoveToCOM(&pLBH->pPOS[0]);
		///////////////////////////////
		size_t ig1=0, ig2=0;
		for(size_t i=0;i<pL->size();i++)
		  {
		    if(pL->R(i)<1.0)
		      {
			ig1++;
			
		      }
		    if(pL->R(i)<2.0)
		      {
			ig2++;
		      }
		  }
		///////////////////////////////
		cout.precision(10);
		cout<<std::setw(12)<<std::fixed<<GetISnap(pL->m_fname)<<" "<<ig1<<" "<<ig2<<" "<<pLBH->GetBHMass()<<endl;		
		///////////////////////////////
		delete pL;
		delete pLBH;
		}


	return 0;
	}




