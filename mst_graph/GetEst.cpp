#include "GetEst.h"
#include <utility>
#include <cstring>
#include <cmath>
#include "kdtree2.hpp"
#include "functions.h"
#include "readers.h"

#include "cats.h"


#define ND_GROUPS 6
GetEst::GetEst()
	{
	dim=ND_GROUPS;
	}

GetEst::~GetEst(void)
	{

	}

template<class T>
void GetEst::fill_ngb_group(int i,boost::numeric::ublas::matrix<T> &ngbGroup,int count, eDOSORT which_est)
	{
	int j=0;
	for(int ingb=0;ingb<count;ingb++)
		{	
		j=Part[i].getNGB(ingb, which_est);
		for(size_t idx=0;idx<ngbGroup.size1();idx++)
			{
			ngbGroup(idx,ingb)=Part[j].Pos[idx];
			}
		}
	}
template<class T>
void GetEst::fill_ngb_group_Ap(int i,boost::numeric::ublas::matrix<T> &ngbGroup,int count, eDOSORT which_est)
	{
	int j=0;
	for(int ingb=0;ingb<count;ingb++)
		{
		float Ap[4];
		j=Part[i].getNGB(ingb, which_est);
		/*for(size_t idx=0;idx<ngbGroup.size1();idx++)
			{
			ngbGroup(idx,ingb)=Part[j].Pos[idx];
			}
			*/
		GetAP(&Part[j].Pos[0], &Ap[0]);
		for(size_t idx=0;idx<ngbGroup.size1();idx++)
			ngbGroup(idx,ingb)=Ap[idx];
		}
	}
void GetEst::run_byKdTree()
	{

	int i=0;
	kdtree2_array realdata;

	realdata.resize(boost::extents[All.NumPart][dim]);
	for (i=0; i<All.NumPart; i++) {
		for (int j=0; j<dim; j++) 
			realdata[i][j] = Part[i].Pos[j];
		}
	std::vector<MyFloat> qv;
	int iq;
	kdtree2_result_vector ngblist, ngblistbrut;
	std::cout << "Build the KD-Tree:.." ;
	kdtree2*  tree = new kdtree2(realdata);
	tree->sort_results=true;
	boost::numeric::ublas::matrix<float> ngbGroup(ND_GROUPS,NUM_NGB);	
	CKernel<float> *pKernel=new CEpanechikov<float>(ND_GROUPS);// Init Kernel for calculations
	int j=0;
	qv.resize(dim);
	for(i=0;i<All.NumPart;i++)
		{
		iq=i;

		qv.assign(&Part[iq].Pos[0],&Part[iq].Pos[0]+dim);

		Part[i].pNGB=new int[NUM_NGB];
		tree->n_nearest( qv , NUM_NGB,ngblist);		
		int ingb=0;
		BOOST_FOREACH( kdtree2_result ngbNum, ngblist )
			{	
			j=ngbNum.idx;
			Part[i].pNGB[ingb++]=j;
			cout<<j<<endl;
			}
		//fill_ngb_group<float>(i,ngbGroup);
		//Metrics::MahalDistance<float> Mahal(ngbGroup);
		float  dist=0.0;
		float hsml2=(float)ngblist[NUM_NGB-1].dis;//Mahal.doDist(&Part[i].Pos[0],&Part[ngblist[NUM_NGB-1].idx].Pos[0]);
		float hsml=std::sqrt(hsml2);
		float u=0.0f, mahalDet=0.0f;
		Part[i].Est=0.0;
		BOOST_FOREACH( kdtree2_result ngbNum, ngblist )
			{	
			j=ngbNum.idx;
			dist=(float)ngbNum.dis;//Mahal.doDist(&Part[i].Pos[0],&Part[j].Pos[0]);
			Part[i].Est+=(pKernel->W(std::sqrt(dist/hsml2)));

			mahalDet+=(float)ngbNum.metric->getDet();
			}
		float hd=1.0;
		for(int idim=0;idim<ND_GROUPS;idim++)
			{
			hd*=hsml;
			}
		Part[i].Est/=(hd*std::sqrt(mahalDet/(float)NUM_NGB));

		//				if(i==5)dumpNgbTofile("c:/arm2arm/DATA/ngb.txt", &Part[i].pNGB[0], NUM_NGB);
		//cout<<i<<" " <<Part[i].Est<<endl;
		}
	/////////////////////////////////////////
	DumpVectorEst("c:/arm2arm/DATA/smooth64.est");
	DumpVectorEst("c:/arm2arm/DATA/smooth64.rho",1);
	/////////////////////////////////////////
	}


// // should be GetAp<CCoord, float>(vec, A)
// A returns  2 elements
template<class T >
void GetAP(T *v, T *A)
	{
	T R=std::sqrt(v[0]*v[0]+v[1]*v[1]);
	T cos_phi=v[0]/R , sin_phi=v[1]/R;
	//Ar=Ax*cos_phi +Ay*sin_phi;
	//Ap=-Ax*sin_phi+Ay*cos_phi;
	A[0]=R;
	A[1]=v[3]*cos_phi +v[4]*sin_phi;
	A[2]=-v[3]*sin_phi+v[4]*cos_phi;
	}

////////////////////////////////////////////////////
template <class T>
void GetEst::DumpVector(std::string fname, std::vector<T> &vec)
	{
	std::ofstream vec_dump(fname.c_str(), std::ios::binary);
	size_t np=vec.size();
	vec_dump.write((char*)&np,sizeof(size_t));
	for(size_t i=0;i<np;i++)
		vec_dump.write((char*)&vec[i],sizeof(T));
	vec_dump.close();
	}
template <class T>
void GetEst::LoadDumpVector(std::string fname, std::vector<T> &vec)
	{
	std::ifstream vec_dump(fname.c_str(), std::ios::binary);
	size_t np=0;
	vec_dump.read((char*)&np, sizeof(np));
	vec.resize(np);	
	for(size_t i=0;i<np;i++)
		vec_dump.read((char*)&vec[i],sizeof(T));
	vec_dump.close();
	}

////////////////////////////////////////////////////
void GetEst::DumpVectorEst(string filename, int which)
	{
	std::ofstream of(filename.c_str(),std::ios::out | std::ios::binary);
	if(!of.good())assert("Error");
	unsigned int i=0, np=All.NumPart;
	of.write((char*)&np,sizeof(uint));
	if(which==0)
		for(i=0;i<np;i++)
			of.write((char*)(&Part[i].Est),sizeof(float));
	else
		for(i=0;i<np;i++)
			of.write((char*)(&Part[i].Rho),sizeof(float));

	of.close();
	}


void GetEst::run()
	{
	int i=0;
	std::vector<float> Rho(All.NumPart),Est(All.NumPart);
	CKernel<float> *pKernel=new CEpanechikov<float>(dim);// Init Kernel for calculations
	eDOSORT which_est=BY_EST;
	int j=0, ihsml=0;
	float u=0.0, hsml=0.0;
	int ngbmax=10;
	int ngb_for_mahal=min(40,All.DesNumNgb); 
	dim=3;
	for(i=0;i<All.NumPart;i++)
		{		
		boost::numeric::ublas::matrix<float> ngbGroup(dim,ngb_for_mahal);	
		fill_ngb_group_Ap<float>(i,ngbGroup,ngb_for_mahal, which_est);
		Metrics::MahalDistance<float> Mahal(ngbGroup);
		Est[i]=0.0;
		Rho[i]=0.0;
		ihsml=Part[i].getNGB(ngbmax,which_est);
		hsml=Mahal.doDist(&Part[i].Pos[0],&Part[ihsml].Pos[0]);
		for( int ingb=0;ingb<ngbmax; ingb++)
			{	
			 j=Part[i].getNGB(ingb,which_est);
			 u=Mahal.doDist(&Part[i].Pos[0],&Part[j].Pos[0]);
			 u/=hsml;
			 Rho[i] +=  pKernel->W(u);
			 Est[i] += (pKernel->W(u)*Part[j].Est);
			 //cout<<j<<" "<<Rho[i]<<" u: "<<u<<" "<<hsml<<endl;
			}
		Est[i]=Rho[i];
		Rho[i]=Part[i].Rho;
		}
	/////////////////////////////////////////
	DumpVector<float>("c:/arm2arm/DATA/smooth64.est", Est);
	DumpVector<float>("c:/arm2arm/DATA/smooth64.rho",Rho);
	/////////////////////////////////////////
	
	
	}

void GetEst::run_graph_Ap()
	{
	Graph graphFOF;
	int i=0;
	std::vector<float> Rho(All.NumPart),Est(All.NumPart), vec(All.NumPart);
	eDOSORT which_est=BY_EST;
	int j=0, ihsml=0;
	float u=0.0, hsml=0.0;
	int ngbmax=All.DesNumNgbA;
	int ngb_for_mahal=All.DesNumNgb;//min(40,All.DesNumNgb);
//	int jj;
	float uj=0.0,dist,wl,Thr;
	dim=3;
	ARR.resize(All.NumPart);
	float mi=1e10,ma=-1e10,meanEst=0.0;
	boost::numeric::ublas::matrix<float> ngbGroup(dim,All.NumPart);	
	for(i=0;i<All.NumPart;i++)
		{
		GetAP(&Part[i].Pos[0], &ARR[i].Pos[0]);
		for(size_t idx=0;idx<ngbGroup.size1();idx++)
			//ngbGroup(idx,i)=Part[i].Pos[idx];//ARR[i].Pos[idx];
			ngbGroup(idx,i)=ARR[i].Pos[idx];
		wl=log10(Part[i].get());
		meanEst+=wl;
		mi=min(wl, mi);
		ma=max(wl, ma);
		}
	meanEst/=(float)All.NumPart;
	

	//boost::numeric::ublas::matrix<float> ngbGroup(dim,ngb_for_mahal);	
		//fill_ngb_group_Ap<float>(i,ngbGroup,ngb_for_mahal, which_est);
	Metrics::MahalDistance<float> Mahal(ngbGroup);
	int jj=0;	
	for(i=0;i<All.NumPart;i++)
		{
		
		//boost::numeric::ublas::matrix<float> ngbGroup(dim,ngb_for_mahal);	
		//fill_ngb_group_Ap<float>(i,ngbGroup,ngb_for_mahal, which_est);
		//Metrics::MahalDistance<float> Mahal(ngbGroup);
		
		//j=Part[i].getNGB(0,which_est);
		//cout<<i<<" "<<j<<endl;
		//cout<<Part[i].Pos[0]<<" "<<Part[i].Pos[1]<<" "<<Part[i].Pos[2]<<endl;
		//u=Mahal.doDist(&Part[i].Pos[0],&Part[j].Pos[0]);
		//u=Mahal.doEDist(&Part[i].Pos[0],&Part[j].Pos[0]);
		//jj=j;
		for( int ingb=0;ingb<ngbmax; ingb++)
			{	
			 j=Part[i].getNGB(ingb,BY_AEST);			 
			// uj = std::sqrt(Mahal.doDist(&ARR[i].Pos[0],&ARR[j].Pos[0]));	
			 //uj = Mahal.doEDist(&Part[i].Pos[0],&Part[j].Pos[0]);	
			 ///if(uj<u && i!=j)
			//	 {
			//	  u=uj;
			//	  jj=j;
			//	 }
			 uj=Part[j].get();
			 //if(ARR[i].Pos[0])
			//	add_edge(i,j, uj,graphFOF);
			 dist=abs(ARR[i].Pos[0]-ARR[j].Pos[0]);
			 Thr=(log10(uj)-meanEst)/(ma-meanEst);
			if( dist<0.1 && Thr>0.2)
			add_edge(i,j, dist,graphFOF);
			}
//		vec[i]=Part[i].get();
//		if(abs(ARR[i].Pos[0]-ARR[jj].Pos[0])<0.1 && vec[i])
//			add_edge(i,jj, u,graphFOF);
		cout<<i<<"\r";
		}
	cout<<endl;
///////////////////////////////////////////////////////////////////
	cout<<"Num Forest edges: "<<num_edges(graphFOF)<<endl;
	std::vector<int> component(num_vertices(graphFOF));
	int num = connected_components(graphFOF, &component[0]);
	cout<<"Number of components:"<<num<<endl;
	cout<<"======================="<<endl;
	/////////////// BUILD The MST tree
	weight_map_type weight = get(edge_weight, graphFOF);
	std::vector < Edge > spanning_tree;

	std::cout << "Building MST:." ;	
	kruskal_minimum_spanning_tree(graphFOF, std::back_inserter(spanning_tree));
	std::cout <<"MST num edges="<< spanning_tree.size()<<" .. done " << std::endl;
////////////////////////////////////////////////////////////////////
//start to CUT
/////////////////////////////	

	DumpMstCat(graphFOF, spanning_tree);

//////////////////////////////
/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////
	DumpVector<float>("c:/arm2arm/DATA/smooth64.est", Est);
	DumpVector<float>("c:/arm2arm/DATA/smooth64.rho",Rho);
	DumpVector<float>("c:/arm2arm/DATA/mahal.u",vec);
	DumpVector<tARR>("c:/arm2arm/DATA/Arr.u",ARR);
	/////////////////////////////////////////
	}

template <class T>
T distance(T *P1,T *P2, int dim=3)		 
	{
		T dist=0.0;
		for(int i=0;i<dim;i++)
			dist+=(P1[i]-P2[i])*(P1[i]-P2[i]);
		dist=std::sqrt(dist);
		return dist;
	}
void GetEst::run_graph()
	{
	Graph graphFOF;
	int i=0;
	std::vector<float> Rho(All.NumPart),Est(All.NumPart), vec(All.NumPart);
	std::vector<uint> isdone(All.NumPart);
	eDOSORT which_est=BY_EST;
	int j=0, ihsml=0;
	float u=0.0, hsml=0.0, esti,estj=0.0, dist;
	int ngbmax=10;//All.DesNumNgbA;	
	cout<<"\nWe are linking particles with "<<All.DesNumNgbA<<endl;
	for(i=0;i<All.NumPart;i++)
		{
		isdone[i]=-1;
		}

	for(int ip=0;ip<All.NumPart;ip++)
		{
		i=isortEst[ip];
		esti=Part[i].get(which_est);
		if(isdone[i]==-1)
		for(int ingb=0;ingb<ngbmax; ingb++)
			{	
			j=Part[i].getNGB(ingb,which_est);
			if(isdone[j]==-1)
				{
				estj=Part[j].get(which_est);
				dist=distance<float>(&Part[i].Pos[0],&Part[j].Pos[0]);		 
				if( estj<esti && dist <0.2)
				 {
				 add_edge(i,j, dist,graphFOF);
				 isdone[j]=i;
				 }
				}
			}

	if(i%10 ==0)
			{
			std::cout<<std::setprecision(2)<<std::fixed;
			std::cout<<ip/double(All.NumPart)*100.0<<"%"<<"\r";
			std::cout<<std::setprecision(8)<<std::fixed;
			}
		}
	cout<<endl;
///////////////////////////////////////////////////////////////////
	cout<<"Num Forest edges: "<<num_edges(graphFOF)<<endl;
	std::vector<int> component(num_vertices(graphFOF));
	int num = connected_components(graphFOF, &component[0]);
	cout<<"Number of components:"<<num<<endl;
	cout<<"======================="<<endl;
	/////////////// BUILD The MST tree
	weight_map_type weight = get(edge_weight, graphFOF);
	std::vector < Edge > spanning_tree;

	std::cout << "Building MST:." ;	
	kruskal_minimum_spanning_tree(graphFOF, std::back_inserter(spanning_tree));
	std::cout <<"MST num edges="<< spanning_tree.size()<<" .. done " << std::endl;
////////////////////////////////////////////////////////////////////
//start to CUT
/////////////////////////////	

	DumpMstCat(graphFOF, spanning_tree);

//////////////////////////////
/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////
	//DumpVector<float>("c:/arm2arm/DATA/smooth64.est", Est);
	//DumpVector<float>("c:/arm2arm/DATA/smooth64.rho",Rho);		
	/////////////////////////////////////////
	}
void GetEst::DumpMstCat(Graph &graphFOF,std::vector < Edge > &spanning_tree)
	{
	
		Graph graphMST;
	for (std::vector < Edge >::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) {
			add_edge(source(*ei, graphFOF),target(*ei, graphFOF),get(edge_weight, graphFOF, *ei),graphMST);
		}
	uint numvert=num_vertices(graphMST);
	std::vector<int> component(numvert);
	int num = connected_components(graphMST, &component[0]);
	cout<<"==============\nTotal Vertexes="<<numvert<<endl;
	cout << "Total number of connected components: " <<std::setprecision(5)<< num <<endl;
	//// make catalogues 	
	CCatalog c;
	std::vector<int>::size_type i;
	c.cats.resize(num);
	for (i = 0; i != component.size(); ++i)
		{
		c.cats[component[i]].id.push_back(i);
		}
	for (i = 0; i != component.size(); ++i)
		c.cats[component[i]].m_Np=c.cats[component[i]].id.size();

	c.sort();
	c.SaveCats("c:/arm2arm/DATA/test.idx");
	}


#include <ANN/ANN.h>					// ANN declarations
void printPt(std::ostream &out, ANNpoint p, int dim=2)			// print point
{
	out << "(" << p[0];
	for (int i = 1; i < dim; i++) {
		out << ", " << p[i];
	}
	out << ")\n";
}

/*	std::ofstream out_dump_file("c:/arm2arm/DATA/kdtree.dmp");
	kdTree->Dump(ANNtrue, out_dump_file);
	out_dump_file.close();*/

template<class T>
void fill_ngb_byID(int *id,boost::numeric::ublas::matrix<T> &ngbGroup)
	{
	int j=0;
	for(size_t ingb=0;ingb<ngbGroup.size2();ingb++)
		{	
		j=id[ingb];
		for(size_t idx=0;idx<ngbGroup.size1();idx++)
			{
			ngbGroup(idx,ingb)=Part[j].Pos[idx];
			}
		}
	
	}
class CCompare{
public:
	CCompare(std::vector<double> &vec):myvec(vec){};
	const	bool operator() ( int& lhs, int& rhs) {
			return myvec[lhs] < myvec[rhs];
		}
private:
		std::vector<double> &myvec;
	};
void GetEst::AllocateANNTreeStructures(int dim,int nPts, int kd)
	{
	queryPt = annAllocPt(dim);					// allocate query point
	dataPts = annAllocPts(nPts, dim);			// allocate data points
	nnIdx = new ANNidx[kd];						// allocate near neigh indices
	dists = new ANNdist[kd];					// allocate near neighbor dists
	};
void GetEst::FillANNData(int dim)
	{
	for(int i=0;i<All.NumPart;i++)
		{
		for(int j=0;j<dim;j++)
			{
			dataPts[i][j]=Part[i].Pos[j];//rand()/(double)RAND_MAX;
			}
		}
	};
void GetEst::run_ANN(std::vector<float> &annRho,int dim,//Dimensions
	int kd,//Estimation Neighbours
	int bucketsize,//Tuning parameter for ANN speed
	bool verbose
	)
	{
	double uj;
	int			nPts=All.NumPart;					// actual number of data points
	

	AllocateANNTreeStructures(dim, nPts, kd);
	FillANNData(dim);
	kdTree = new ANNkd_tree(					// build search structure
		dataPts,					// the data points
		nPts,						// number of points
		dim,						// dimension of space
		bucketsize,	                // bucket size
		ANN_KD_SUGGEST              //  Splitting rule
		//ANN_KD_STD, standard kd-splitting rule 
		//ANN_KD_MIDPT, midpoint split 
		//ANN_KD_FAIR, fair-split 
		//ANN_KD_SL_MIDPT, sliding midpoint split 
		//ANN_KD_SL_FAIR, sliding fair-split 
		//ANN_KD_SUGGEST the authors' suggestion for best
		
		);

///////// START to search NGBs ////
	annRho.resize(nPts);
  	std::vector<double> dist(kd);
	CKernel<double> *pKernel=new CEpanechikov<double>(dim);// Init Kernel for calculations
	boost::numeric::ublas::matrix<double> ngbGroup(dim,kd);	
	Metrics::MahalDistance<double> Mahal;
	std::vector<int> indns, ind;
	indns.resize(kd);
	for(int i=0;i<kd;i++)
		indns[i]=i;

	double hsml;
	if(dim==6)//get 6D NGB
		{
		Klink=20;
		Mngb.resize(nPts);
		 for(size_t i=0;i<Mngb.size();i++)
			 Mngb[i].p.resize(Klink);
		}
	for(int i=0;i<All.NumPart;i++)
		{
		kdTree->annkSearch(						// search
			&dataPts[i][0],						// query point
			kd,								// number of near neighbors
			nnIdx,							// nearest neighbors (returned)
			dists							// distance (returned)
			//eps // error bound
			);								
		fill_ngb_byID<double>(nnIdx,ngbGroup);
		Mahal.Init(ngbGroup);
		
		for (int j = 0; j < kd; j++) {
			uj=std::sqrt(Mahal.doDist(&Part[i].Pos[0],&Part[nnIdx[j]].Pos[0]));
			dist[j]=uj;
			}
		if(dim==6)
			{
			CCompare functor_idcomp(dist);
			ind=indns;
			std::sort(ind.begin(),ind.end(), functor_idcomp);
			for(int kl=0;kl<Klink;kl++)
				Mngb[i].p[kl]=std::make_pair<int, float>(nnIdx[ind[kl]],(float)dist[ind[kl]]);
			}
			std::sort(dist.begin(),dist.end());

		hsml=dist[kd-1];
		for (int j = 0; j < kd; j++) {
			uj=dist[j];
			uj/=hsml;
			annRho[i]+=(float)pKernel->W(uj);
			}
		double hd=1.0;
		for(int idim=0;idim<dim;idim++)
			{
			hd*=hsml;
			}

		annRho[i]/=float(hd*std::sqrt(Mahal.getDet()));
	
		
			if(verbose && (i%100==0))
				cout<<i<<"\r";
		}
//////////////// DONE analysis  ////////////////////////////
	delete pKernel;
	DeallocateANN();
	}
void GetEst::Run_SPHEst()
	{
	vector<float> annRho;
	vector<float> annEst;
	vector<float> annSmooth;
	
		if(false){ 
		scoped_timer timethis("#GetEst::Run_SPHEst():> 3D density: run_ANN\t"); 

		run_ANN(annRho, 3,//Dimensions
			40,//Estimation Neighbours
			3-1,//Tuning parameter for ANN speed
			true
			);	
		}else
			LoadDumpVector<float>("c:/arm2arm/DATA/smooth64.rho",annRho);

	if(true)
		{ 
		scoped_timer timethis("#GetEst::Run_SPHEst():> 6D density: run_ANN\t"); 

		run_ANN(annEst, 6,//Dimensions
			64,//Estimation Neighbours
			6-1,//Tuning parameter for ANN speed
			true
			);	
		DumpNgb("c:/arm2arm/DATA/Mngb.est");
		}else{
			LoadDumpVector<float>("c:/arm2arm/DATA/smooth64.est",annEst);
			//LoadNGB("c:/arm2arm/DATA/Mngb.est");
			}

		if(false){ 
		scoped_timer timethis("#GetEst::Run_SPHEst:> Smoothing Est: SmoothByANN\t"); 
		
		
		
		SmoothByANN(annRho,//Based on what to smooth
			annEst,//What to smooth
			annSmooth,//result
			5//Smoothing Neighbours in 3D
			);	
		DumpVector<float>("c:/arm2arm/DATA/smooth64.sme",annSmooth);
		}
		

		DumpVector<float>("c:/arm2arm/DATA/smooth64.rho",annRho);
		DumpVector<float>("c:/arm2arm/DATA/smooth64.est",annEst);
		

	}

void GetEst::LoadNgb(string fname)
	{
	std::ifstream vec_dump(fname.c_str(), std::ios::binary);
	size_t np, kl;
	
	vec_dump.read((char*)&np,sizeof(size_t));
	vec_dump.read((char*)&kl,sizeof(size_t));
	
	Mngb.resize(np);
	
	for(size_t i=0;i<np;i++)
		{
		Mngb[i].p.resize(kl);
		for(size_t ki=0;ki<kl;ki++)
			{
			vec_dump.read((char*)&Mngb[i].p[ki].first,sizeof(int));
			vec_dump.read((char*)&Mngb[i].p[ki].second,sizeof(float));
			}
		}
		vec_dump.close();
	}

void GetEst::DumpNgb(std::string fname)
	{
	std::ofstream vec_dump(fname.c_str(), std::ios::binary);
	size_t np=Mngb.size(), kl=Mngb[0].p.size();
	vec_dump.write((char*)&np,sizeof(size_t));
	vec_dump.write((char*)&kl,sizeof(size_t));
	for(size_t i=0;i<np;i++)
		for(size_t ki=0;ki<kl;ki++)
			{
			vec_dump.write((char*)&Mngb[i].p[ki].first,sizeof(int));
			vec_dump.write((char*)&Mngb[i].p[ki].second,sizeof(float));
			}
	vec_dump.close();
	}
void GetEst::makeMSTree()
	{
	Graph graphFOF;
	float dist;
	int j;
////////////////////////////////////////////////////////////////////////
		for(size_t i=0;i<Mngb.size();i++)
			{
			for(size_t kl=0;kl<Mngb[i].p.size();kl++)
				{
				dist=Mngb[i].p[kl].second;
				j=Mngb[i].p[kl].first;
				if(dist<1.1)
					add_edge(i,j, dist,graphFOF);
				}
			}
////////////////////////////////////////////////////////////////////////

	cout<<"Num Forest edges: "<<num_edges(graphFOF)<<endl;
	std::vector<int> component(num_vertices(graphFOF));
	int num = connected_components(graphFOF, &component[0]);
	cout<<"Number of components:"<<num<<endl;
	cout<<"======================="<<endl;
	/////////////// BUILD The MST tree
	weight_map_type weight = get(edge_weight, graphFOF);
	std::vector < Edge > spanning_tree;

	std::cout << "Building MST:." ;	
	kruskal_minimum_spanning_tree(graphFOF, std::back_inserter(spanning_tree));
	std::cout <<"MST num edges="<< spanning_tree.size()<<" .. done " << std::endl;
	////////////////////////////////////////////////////////////////////
	//start to CUT
	/////////////////////////////	

	DumpMstCat(graphFOF, spanning_tree);

	
	}
void GetEst::SmoothByANN(std::vector<float> &annRho, std::vector<float> &annEst,//What to smooth
						 std::vector<float> &annSmooth,//result
						 int nsph,//Smoothing Neighbours in 3D
						 bool verbose
						 )
	{
	int			nPts=All.NumPart;// actual number of data points
	int dim=3;//dims
	annSmooth.resize(nPts);
	AllocateANNTreeStructures(dim, nPts, nsph);
	FillANNData(dim);
	kdTree = new ANNkd_tree(					// build search structure
		dataPts,					// the data points
		nPts,						// number of points
		dim//,						// dimension of space
		//bucketsize,	                // bucket size
		//ANN_KD_SUGGEST              //  Splitting rule
		//ANN_KD_STD, standard kd-splitting rule 
		//ANN_KD_MIDPT, midpoint split 
		//ANN_KD_FAIR, fair-split 
		//ANN_KD_SL_MIDPT, sliding midpoint split 
		//ANN_KD_SL_FAIR, sliding fair-split 
		//ANN_KD_SUGGEST the authors' suggestion for best

		);

	CKernel<double> *pKernel=new CEpanechikov<double>(dim);// Init Kernel for calculations	
		
	double hsml=0.0,uj=0.0;
	double hd;
	double sumR=0.0;
	double sumRho=0.0;
	double PI_COEF=M_PI*4.0/3.0;

	for(int i=0;i<All.NumPart;i++)
		{
		kdTree->annkSearch(						// search
			&dataPts[i][0],						// query point
			nsph,								// number of near neighbors
			nnIdx,							// nearest neighbors (returned)
			dists							// distance (returned)			
			);								
	
		hsml=dists[nsph-1];
		sumR=0.0;
		sumRho=0.0;
		for (int j = 0; j < nsph; j++) {
			uj=dists[j];
			uj/=hsml;
			sumR+=(float)pKernel->W(uj)*annEst[nnIdx[j]];
			sumRho+=(float)pKernel->W(uj);
			}
		hd=hsml*hsml*hsml;		
		annSmooth[i]=(float)(sumR/(hd*PI_COEF*sumRho));
		if(verbose && (i%100==0))
			cout<<i<<"\r";
		}
	//////////////// DONE analysis  ////////////////////////////
	delete pKernel;
	DeallocateANN();
	};	