// profile.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>
#include <numeric>
#include <valarray>
#include <functional>
#include "loader.h"
#include "options.h"
#include "utils.h"
#include "Logger.h"
#include "MSTree.h"
#include "PCA.h"
#include "OAFHelper.h"
#include "Sigma.h"
using std::string;

//////////////////////
//Result format R1, Am1, R2,  Am2, R3, Am3
std::vector<float> getAmMaxMeanR(TLogData &data){ 
	std::vector<float> result;

	//find the maximum id;
	vector<float> R=data[0];
	unsigned int Rshift=std::distance(
		std::find_if(R.begin(), R.end(), std::bind2nd( std::greater_equal<float>(), 10 )), 
		R.end()
		);
	vector<float> Am=data[1];
	vector<float>::iterator itmax=std::max_element(Am.begin(), Am.end());
	if(itmax==Am.end())
		{ cout<<"Error to find the maximum"<<endl;return result;}
	float Am1=(0.5f*(*itmax));
	vector<float>::iterator it2=std::find_if(itmax, Am.end(), std::bind2nd( std::less_equal<float>(), Am1 )  );
	if(it2==Am.end()){ cout<<"Error to find the Half of maximum"<<endl;return result;}
	result.push_back(R[std::distance(Am.begin(),itmax)]);
	result.push_back(*itmax);
	result.push_back(R[std::distance(Am.begin(),it2)]);
	result.push_back( (*it2));
	float val= *it2;
	while( *(++it2) < val)// && (it2+1)!=(Am.end()-1))
		{ val= (*it2);};
	if(it2==Am.end()){ cout<<"Error to find the next Minimum of maximum"<<endl;return result;}
	result.push_back(R[std::distance(Am.begin(),it2)]);
	result.push_back( (*it2));

	return result;
	};

////////////////////////////////////////////

////////////////////////////////////////////

int main(int argc, char* argv[])
	{
	scoped_timer timemme("Main program Profiler :.....");
	COptions opt(argc, argv);
	///////////////////////// We will move to another directory
	typedef boost::ptr_vector<COAFHelper> TOAFHelper_vec;
	TOAFHelper_vec vecOAFABC;
	if(opt.m_OAF)
		{
		vecOAFABC.push_back(new COAFHelper("bulge"));
		vecOAFABC.push_back(new COAFHelper("disk"));
		vecOAFABC.push_back(new COAFHelper("halo"));
		for(size_t i=0;i<opt.m_IDlistvec.size();i++)
			vecOAFABC[i].SetID(opt.m_IDlistvec[i]);
		}
/////////////////////////
	if(opt.is_bad())
		return EXIT_FAILURE;
	unsigned int ParticleType=opt.m_type;
	typedef std::vector<CLoader*> TLoader;
	TLoader Lvec;
	std::vector<float> rvec;
	const int NFields=5+5+4+1;
	CLogger log(opt.m_file_out,opt.m_updatelog, NFields);
	CLogger logAmRR("AmRad.log", opt.m_updatelog);
	logAmRR.SetComment("Am R1 R2 R3");
	CLoggerAm logAm("Amf.log",opt.m_updatelog);
	logAm.SetComment("R= 0.1-20.0");

	log.SetComment(" snap  r1kpc(0-4) r2kpc(0-4) ");
	rvec.push_back(1.0f);
	rvec.push_back(2.0f);

	int ir=0;
	for(size_t i=0;i<opt.m_snapshotList.size();i++)
		{
		unsigned int isnap=GetISnap(opt.m_snapshotList[i]);  
		if(!is_file_exist(opt.m_snapshotList[i])){cout<<"skipping...isnap="<<isnap<<endl;continue;}

		if(log.is_done(isnap))continue;
		CLoader *pL=new CLoader(opt.m_snapshotList[i],ParticleType, true);
                 cout<<"We got Np:"<<pL->size()<<endl;
		if(opt.m_OAF)
			{
			opt.m_OAF=pL->ReadID();
			
			}
		if(true)
			{
			///////////////////////////////
			// Do profiling here
			std::vector<float> fVec(NFields), x, y,z, R;
			vector<unsigned int > idxR(pL->m_nelem-1,0);	
			///////////////////////////////////
			//cout<<"COM: "<<pL->m_COM[0]<<" "<<pL->m_COM[1]<<" "<<pL->m_COM[2]<<endl;
			CRange range[3];
			for(unsigned int i=0, ig=0, ist=0;i<pL->m_nelem-1;i++)// -1 to exclude BH particle
				if( pL->pType[i] == 4 )
					{
					x.push_back(pL->pPOS[i*3]);
					y.push_back(pL->pPOS[i*3+1]);
					z.push_back(pL->pPOS[i*3+2]);	
					range[0].getbound(x[ist]);
					range[1].getbound(y[ist]);
					range[2].getbound(z[ist]);
					ist++;
					}
				range[0].print("# Range for x coord");
				range[1].print("# Range for y coord");
				range[2].print("# Range for z coord");
				cout<<"# got a :"<<x.size()<<" particles for determining the center"<<endl;
				
				CMSTree mst_tree(x,y,z, opt.m_eps, opt.m_min_npart, opt.m_NGB);
				//mst_tree.m_MSTCatalog[0].wcom[0];
			    float com[]={ 19289.00659f, 26536.94112f, 24105.424500f};
				com[0]=(float)mst_tree.m_MSTCatalog[0].wcom[0];
				com[1]=(float)mst_tree.m_MSTCatalog[0].wcom[1];
				com[2]=(float)mst_tree.m_MSTCatalog[0].wcom[2];
				
				pL->MoveToCOM(&com[0]);
				std::transform( x.begin(), x.end(), x.begin(),std::bind2nd( std::minus<float>(), com[0]) );
				std::transform( y.begin(), y.end(), y.begin(),std::bind2nd( std::minus<float>(), com[1]) );
				std::transform( z.begin(), z.end(), z.begin(),std::bind2nd( std::minus<float>(), com[2]) );
			
				R.clear();
				for(unsigned int i=0, ig=0, ist=0;i<pL->m_nelem-1;i++)// -1 to exclude BH particle
					{
					float rr=pL->R(i);//

					if(pL->pType[i] == 0)ig++;
					if(pL->pType[i] == 4)
						{
						R.push_back(sqrt(x[ist]*x[ist]+y[ist]*y[ist]+z[ist]*z[ist]));
						ist++;
						}

					if(opt.m_OAF)
					if(pL->pType[i]==0 || pL->pType[i]==4)
						{
							{
							for(size_t iv=0;iv<vecOAFABC.size();iv++)
								{
								float mass, sfr, hsml, rho, u, zm[2];
								unsigned char type=pL->pType[i];
								
								sfr=pL->pSFR[ig-1];
								hsml=pL->pHSML[ig-1];
								rho=pL->pRHO[ig-1];
								
								mass=pL->pMASS[ig-1+ist-1-1];
								u=pL->pU[ig-1];
								if(type==0)zm[0]=pL->pZ[ig-1];
								else
									zm[1]=pL->pZ[ig-1+ist-1-1];
							
//3Pos, 3Vel, MASS, HSML, RHO, U, Z, SFR							
								if(vecOAFABC[iv].insert(
									pL->pID[i],isnap,type, 
									&pL->pPOS[i*3], 
									&pL->pVEL[i*3],
									mass,hsml, rho, u, &zm[0], sfr)
									)
									break;
								}
							}

						}


					for(ir=0;ir<2;ir++)
						if(rr<rvec[ir]){
							fVec[ir*5+pL->pType[i]]+=1.0;
							if(pL->pType[i]==0){
								fVec[10+ir]+=pL->pSFR[ig];
								}
							}

					}
				///////////////////////////////////
				if(true)	{
					CSigma<double> sigma( "234", &pL->pType[0],&pL->pPOS[0], &pL->pVEL[0],pL->size());
					sigma.m_fname="sigma_"+boost::lexical_cast<std::string>(isnap)+string("_4.txt");
					}
///////////////////////////////////////////////////////
				TLogData result=GetAB<float>(R,x,y, 8, 0.1f);//opt.m_Rmax);
				std::vector<float> AmRad=getAmMaxMeanR(result);//format R1, Am1, R2,  Am2, R3, Am3

				/////////////////////////////////		     
				CPCA pca;
				double phi=pca.GetCovarMatrix(mst_tree, 0);
				fVec[14]=(float)phi;
				mst_tree.dump(0);
				logAm.insert(isnap, result);		
				logAmRR.insert(isnap, AmRad);
				log.insert(isnap, fVec);
				pca.dump();
				///////////////////////////////
			}
		delete pL;
		//Lvec.push_back(pL);
		}

	/// Lets free the memory
	//delete_pointers_list< TLoader >(Lvec);

	return 0;
	}




