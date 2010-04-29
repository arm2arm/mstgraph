// profile.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//

#include "loader.h"
#include "options.h"
#include "utils.h"
#include "Logger.h"

int main(int argc, char* argv[])
	{
	COptions opt(argc, argv);
	if(opt.is_bad())
		return EXIT_FAILURE;
	unsigned int ParticleType=opt.m_type;
	typedef std::vector<CLoader*> TLoader;
	TLoader Lvec;
	std::vector<float> rvec;
	CLogger log(opt.m_file_out,opt.m_updatelog);

	log.SetComment(" snap  r1kpc(0-4) r2kpc(0-4) ");
	rvec.push_back(1.0f);
	rvec.push_back(2.0f);

	int ir=0;
	for(size_t i=0;i<opt.m_snapshotList.size();i++)
		{
		CLoader *pL=new CLoader(opt.m_snapshotList[i],ParticleType);
		(*pL).PrintStat();
		///////////////////////////////
		// Do profiling here
		unsigned int isnap=GetISnap(opt.m_snapshotList[i]);
		std::vector<float> fVec(5+5+4);
		///////////////////////////////////
		for(unsigned int i=0, ig=0, ist=0;i<pL->m_nelem-1;i++)// -1 to exclude BH particle
			{
			
			float rr=pL->R(i);
			if(pL->pType[i] == 0)ig++;
			if(pL->pType[i] == 4)ist++;
			for(ir=0;ir<2;ir++)
				if(rr<rvec[ir]){
					fVec[ir*5+pL->pType[i]]+=1.0;
					if(pL->pType[i]==0){
					 fVec[10+ir]+=pL->pSFR[ig];
						}else
							if(pL->pType[i]==4){
								fVec[12+ir]+=pL->pSFR[ist+pL->m_npart[0]-1];
								}

					}

			}
		///////////////////////////////////
		log.insert(isnap, fVec);
		///////////////////////////////
		delete pL;
		//Lvec.push_back(pL);
		}

		/// Lets free the memory
		//delete_pointers_list< TLoader >(Lvec);
			
		return 0;
	}

