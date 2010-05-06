// profile.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//

#include "loader.h"
#include "options.h"
#include "utils.h"
#include "Logger.h"
#include "MSTree.h"
#include "PCA.h"
#include <numeric>
#include <functional>
template<class T> 
std::vector<T> smooth(std::vector<T> &v)
	{

	if(v.size()<5)return v;
	int np=v.size();
	std::vector<T> vsm=v;

	for(int i=2;i<np-2;i++)
		vsm[i]=(v[i-2]+2*v[i-1]+3*v[i]+2*v[i+1]+v[i+2])/9.0f;
	vsm[0]=(vsm[0]+vsm[1]+vsm[2])/3.0f;
	vsm[1]=(vsm[1]+vsm[2]+vsm[3])/3.0f;
	vsm[np-2]=(vsm[np-2]+vsm[np-3]+vsm[np-4])/3.0f;
	vsm[np-1]=(vsm[np-1]+vsm[np-2]+vsm[np-3])/3.0f;
	return vsm;
	}

template<class T>
std::vector< std::vector<T> > GetAB(vector<T> &R,vector<T> &x,vector<T> &y, T profRMAX=10.0, T xStep=0.1)
	{
	std::vector<unsigned int> idxR(R.size(), 0);
	for(size_t ii=0;ii<R.size(); ii++)
		{
		idxR[ii]=ii;
		}

	std::sort(idxR.begin(), idxR.end(), ProxyLess< vector<float> >(R)); //sorting by Radius;

	/////////////////////////////////////////////
	//GetAB
	std::vector< vector<T> > result;
	int fm, ii;
	int iR20 = std::count_if(R.begin(), R.end(),
		std::bind2nd(std::less_equal<T>(), profRMAX));

	int nshell = static_cast<int> (profRMAX / xStep) + 1;
	for (fm = 2; fm < 9; fm += 2) {
		vector<T> al, al0, al2, bl2, fR;
		al.resize(nshell);
		al2.resize(nshell);
		bl2.resize(nshell);
		al0.resize(nshell);
		fR.resize(nshell);
		std::fill(al.begin(), al.end(), (T)0.0);
		std::fill(al0.begin(), al0.end(),(T) 0.0);
		std::fill(al2.begin(), al2.end(), (T)0.0);
		std::fill(bl2.begin(), bl2.end(), (T)0.0);
		std::fill(fR.begin(), fR.end(), (T)0.0);

		int ishell = 0;
		T theta = 0.0;
		T massp = 1.0;


		for (int i = 0; i < iR20; i++) {
			ii=idxR[i];
			ishell = static_cast<int> (R[ii] / xStep);
			theta = atan2(y[ii], x[ii]);
			al2[ishell] += massp * cos(fm * theta);
			bl2[ishell] += massp * sin(fm * theta);
			al0[ishell] += massp;
			}

		for (int i = 0; i < nshell; i++) {
			al[i] = sqrt(al2[i] * al2[i] +
				bl2[i] * bl2[i]);
			if (al0[i] > 0.0f)
				al[i] /= al0[i];
			fR[i] = xStep*i;
			}
		double total = std::accumulate(al.begin(), al.end(), 0.0);
		//string fmFile = fname + string("_Afm") + ToString(fm) + '_' + ToString(TYPE);
		if(fm == 2)	result.push_back(fR);
		result.push_back(smooth<float>(al));
		//result.push_back(al);
		//result.push_back(al2);
		//result.push_back(bl2);
		}
	return result;
	}

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

////////////////////////
void GetOmegaBar(vector<float> &x, vector<float> &y,vector<float>  &z)
	{
	
	
	}
///////////////////////
int main(int argc, char* argv[])
	{
	COptions opt(argc, argv);
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
		CLoader *pL=new CLoader(opt.m_snapshotList[i],ParticleType);

		///////////////////////////////
		// Do profiling here
		std::vector<float> fVec(NFields), x, y,z, R;
		vector<unsigned int > idxR(pL->m_nelem-1,0);	
		///////////////////////////////////
		for(unsigned int i=0, ig=0, ist=0;i<pL->m_nelem-1;i++)// -1 to exclude BH particle
			{
			float rr=pL->R(i);
			idxR.push_back(i);
			if(pL->pType[i] == 0)ig++;
			if(pL->pType[i] == 4)ist++;
			if( (pL->pType[i] == 4 /*|| pL->pType[i] == 2*/) ){
				R.push_back(rr);
				x.push_back(pL->pPOS[i*3]);
				y.push_back(pL->pPOS[i*3+1]);
				z.push_back(pL->pPOS[i*3+2]);			      
				}
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
		TLogData result=GetAB<float>(R,x,y);
		std::vector<float> AmRad=getAmMaxMeanR(result);//format R1, Am1, R2,  Am2, R3, Am3
		
		
		/////////////////////////////////
		CMSTree mst_tree(x,y,z, opt.m_eps, opt.m_min_npart);
		CPCA pca;
		double phi=pca.GetCovarMatrix(mst_tree, 0);
		fVec[14]=(float)phi;
		//mst_tree.dump(0);
		logAm.insert(isnap, result);		
		logAmRR.insert(isnap, AmRad);
		log.insert(isnap, fVec);
		///////////////////////////////
		delete pL;
		//Lvec.push_back(pL);
		}

	/// Lets free the memory
	//delete_pointers_list< TLoader >(Lvec);

	return 0;
	}




