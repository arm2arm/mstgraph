//=======================================================================
// Copyright 2009 Arman Khalatyan, OAMP/LAM  
// distributed under GPL license
//=======================================================================
#include <utility>
#include <map>
#include <set>
#include <cstring>
#include "mst_graph.h"
//#include "ranker.h"
#include "coord.h"
#include "MSTGroup.h"
#include "kdtree2.hpp"
#include <math.h>
#include "functions.h"

#include "Render.h"
#include "HOP.h"

/////////////////////////////
#define ND_GROUPS 2
#define LOAD_DEBUG 0
#define SMALL_TEST 0
#define MIN_NGRP 50
#define MAXNGB 64.0 
std::string base;//this is used for outputs
#define RAND_FRACTION 1.2
///////////////////////////////

bool less_coord(CCoord *c1,CCoord *c2)
	{
	return (*c1).w < (*c2).w;	 
	}
CCoord gencoord(void )
	{
	static int counter=0;
	CCoord loc_coord(counter,
		rand()/MyFloat(RAND_MAX),
		rand()/MyFloat(RAND_MAX),
		rand()/MyFloat(RAND_MAX));

	counter++;
	return loc_coord;
	};
void print_data(typeVecData vec)
	{
	std::copy(vec.begin(), vec.end(), std::ostream_iterator<CCoord>(std::cout, " "));
	}


//
// define, for convenience a 2d array of floats. 
//  
typedef boost::multi_array<MyFloat,2> array2dfloat;

//This is the main 
void GetAP(CCoord v1, MyFloat *A)
	{
	MyFloat R=std::sqrt(v1.pos[0]*v1.pos[0]+v1.pos[1]*v1.pos[1]);
	MyFloat cos_phi=v1.pos[0]/R , sin_phi=v1.pos[1]/R;

	//Ar=Ax*cos_phi +Ay*sin_phi;
	//Ap=-Ax*sin_phi+Ay*cos_phi;
	A[0]=v1.vel[0]*cos_phi +v1.vel[1]*sin_phi;
	A[1]=-v1.vel[0]*sin_phi+v1.vel[1]*cos_phi;
	}

MyFloat CoordDist(CCoord v1,CCoord v2)
	{
	MyFloat ret = 
		(v1.pos[0]-v2.pos[0])*(v1.pos[0]-v2.pos[0])+
		(v1.pos[1]-v2.pos[1])*(v1.pos[1]-v2.pos[1])+
		(v1.pos[2]-v2.pos[2])*(v1.pos[2]-v2.pos[2]);
	return sqrt(ret);
	}

///////////////////////////////
#include "readers.h"
void read_gadget(string fname, typeVecData *data)
	{
	All.ONLY_TYPE=4;
	read_ic12(fname.c_str());
	MyFloat COM[]={0.0,0.0,0.0}, totEst=0.0;
	for(int i=0;i<All.NumPart;i++)
		{
		COM[0]+=Part[i].Pos[0]*Part[i].Est;		
		COM[1]+=Part[i].Pos[1]*Part[i].Est;		
		COM[2]+=Part[i].Pos[2]*Part[i].Est;		
		totEst+=Part[i].Est;
		}
	COM[0]/=totEst;
	COM[1]/=totEst;
	COM[2]/=totEst;
	for(int i=0;i<All.NumPart;i++)
		{
		Part[i].Pos[0]-=(float)COM[0];
		Part[i].Pos[1]-=(float)COM[1];
		Part[i].Pos[2]-=(float)COM[2];
		}
	}
/////////////////////////////
void read_ascii(string fname, typeVecData *data)
	{
	(*data).clear();
	std::ifstream fin(fname.c_str());
	double x, y, z;
	unsigned int i=0;
	srand(100);
	while(!fin.eof())
		{
		fin>>x>>y>>z;	
		if(rand()/MyFloat(RAND_MAX) < 0.20)
			(*data).push_back(CCoord(i, x, y, z));
		i++;
		}
	fin.close();

	}
void read_ascii_vec(string fname, vector<MyFloat> *data)
	{
	(*data).clear();
	std::ifstream fin(fname.c_str());
	double x;
	unsigned int i=0;
	srand(100);
	while(!fin.eof())
		{
		fin>>x;	
		if(rand()/MyFloat(RAND_MAX) < 0.20)
			(*data).push_back( x);
		i++;
		}
	fin.close();

	}
void read_ascii_vecn(string fname,string fwname, typeVecData *data, int n=6)
	{
	(*data).clear();
	std::ifstream fin(fname.c_str());
	std::ifstream fwin(fwname.c_str());

	unsigned int i=0;
	srand(100);
	CCoord t(i,0,0,0);
	while(!fin.eof())
		{
		fin>>t.pos[0];
		fin>>t.pos[1];
		fin>>t.pos[2];
		fin>>t.vel[0];
		fin>>t.vel[1];
		fin>>t.vel[2];
		fwin>>t.w;
		t.id=i;
		if(rand()/MyFloat(RAND_MAX) < RAND_FRACTION)
			(*data).push_back( t);
		i++;
		}
	fin.close();
	fwin.close();
	cout<<"Got Np= "<<(*data).size()<<endl;
	}
////////////////////////////////////////////////////
void DumpVectorEst(string filename, int which=0)
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
void LoadVectorEst(string filename, int which=0)
	{
	std::ifstream of(filename.c_str(),std::ios::in | std::ios::binary);
	if(!of.good())assert("Error");
	unsigned int i=0, np=All.NumPart;
	of.read((char*)&np,sizeof(uint));
	if(which==0)
		for(i=0;i<np;i++)
			of.read((char*)(&Part[i].Est),sizeof(float));
	else
		for(i=0;i<np;i++)
			of.read((char*)(&Part[i].Rho),sizeof(float));

	of.close();
	}

void SmoothSph()
	{
	if(LOAD_DEBUG)
		{
		LoadVectorEst("c:/arm2arm/DATA/smooth64.est");
		LoadVectorEst("c:/arm2arm/DATA/smooth64.rho",1);
		}else
		{

		typedef boost::multi_array<MyFloat,2> array2dfloat;

		array2dfloat realdata; 
		uint dim=3,N=All.NumPart;
		realdata.resize(boost::extents[All.NumPart][dim]);
		for (uint i=0; i<N; i++) {
			for (uint j=0; j<dim; j++) 
				realdata[i][j] = Part[i].Pos[j];
			}
		kdtree2_result_vector ngblist, ngblistbrut;
		std::cout << "Build the KD-Tree:.." ;
		kdtree2*  tree = new kdtree2(realdata);
		tree->sort_results = true;
		cout<<".DONE."<< std::endl;
		///////////////////////////////
		uint j=0;
		const MyFloat PI43=MAXNGB/(4.0/3.*M_PI); 
		MyFloat smEst=0.0, smRho=0.0, hsml, val=0.0;;
		std::vector<MyFloat> qv;//query vector
		for(uint i=0;i<N;i++)
			{
			uint iq=i, nn=0;
			qv.clear();
			qv.push_back(Part[iq].Pos[0]);
			qv.push_back(Part[iq].Pos[1]);
			qv.push_back(Part[iq].Pos[2]);
			tree->n_nearest( qv , (uint)MAXNGB,ngblist);
			smEst=0.0;smRho=0.0;
			hsml=ngblist[(uint)MAXNGB-1].dis;
			Part[i].pNGBR=new int[(uint)(MAXNGB)];
			BOOST_FOREACH( kdtree2_result ngbNum, ngblist )
				{
				j=ngbNum.idx;
				Part[i].pNGBR[nn]=j;
				val=Wsph<MyFloat>(ngbNum.dis,hsml);
				smRho+=val;
				nn++;
				}

			BOOST_FOREACH( kdtree2_result ngbNum, ngblist )
				{
				j=ngbNum.idx;
				smEst+=(Part[j].Est/smRho*Wsph<MyFloat>(ngbNum.dis,hsml));
				}

			Part[iq].Rho=(float)smRho;
			Part[iq].Est=(float)smEst;
			if(i%10 ==0)
				{
				std::cout<<std::setprecision(2)<<std::fixed;
				std::cout<<i/double(N)*100.0<<"%"<<"\r";
				std::cout<<std::setprecision(16)<<std::fixed;
				}

			}
		DumpVectorEst("c:/arm2arm/DATA/smooth64.est");
		DumpVectorEst("c:/arm2arm/DATA/smooth64.rho",1);
			}//end for SmoothSph() 
	}
template<class T>
void fill_EST_group(boost::numeric::ublas::matrix<T> &ngbGroup, int nmax=20000)
	{
	int j=0;
	int maxnp=std::min(All.NumPart,nmax);
	ngbGroup.resize(ND_GROUPS, maxnp);
	for(int i=0;i<maxnp;i++)
		{	
		j=isortEst[i];
		 for(int idx=0;idx<ND_GROUPS;idx++)
			 ngbGroup(idx,i)=Part[j].Pos[idx];			
		}
	}
template<class T>
void fill_ngb_group(int i,boost::numeric::ublas::matrix<T> &ngbGroup)
	{
	int j=0;
	for(int ingb=0;ingb<All.DesNumNgb;ingb++)
		{	
		 j=Part[i].getNGB(ingb);
		 for(int idx=0;idx<ND_GROUPS;idx++)
			 ngbGroup(idx,ingb)=Part[j].Pos[idx];			
		}
	}
void MahalanobisFOF()
	{
	using namespace boost;
	typedef adjacency_list < vecS, vecS, undirectedS,
		no_property, property < edge_weight_t, MyFloat > > Graph;
	typedef property_map<Graph, edge_weight_t>::type EdgeWeightMap;
	typedef property_map<Graph, edge_weight_t>::type	weight_map_type;
	typedef graph_traits < Graph >::edge_descriptor Edge;
	typedef graph_traits < Graph >::vertex_descriptor Vertex;
	typedef std::pair<int, int> E;
	cout<<"Starting to build FOF Graph based on adaptive Ngb list"<<endl;
	Graph graphFOF;
	int j=0;
	boost::numeric::ublas::matrix<double> ngbGroup(ND_GROUPS,All.DesNumNgb);	
	std::ofstream fout("c:/arm2arm/DATA/edge_len.txt");
	fill_EST_group<double>(ngbGroup);
	Metrics::MahalDistance<double> Mahal(ngbGroup);
	boost::numeric::ublas::vector<double> p(ND_GROUPS);
	for(int i=0;i<All.NumPart;i++)
		{
	
		fill_ngb_group(i,ngbGroup);
		Metrics::MahalDistance<double> Mahal(ngbGroup);

		for(int ingb=0;ingb<All.DesNumNgb;ingb++)//
			{

			j=Part[i].getNGB(ingb);
			//if(j!=i)
				{
				
				for(int idx=0;idx<ND_GROUPS;idx++)p(idx)=Part[j].Pos[idx];

				add_edge(i,j,Mahal.doDist(p),graphFOF);
				}//else
				//	Mahal.last_dist=0.0;
			//cout<<Mahal.last_dist<<" Est="<<Part[j].get()<<endl;
			}

		for(int idx=0;idx<ND_GROUPS;idx++)p(idx)=Part[i].Pos[idx];
		fout<<i<<" "<<j<<" "<<Mahal.doDist(p)<<endl;
		if(i%10 ==0)
			{
			std::cout<<std::setprecision(2)<<std::fixed;
			std::cout<<i/double(All.NumPart)*100.0<<"%"<<"\r";
			std::cout<<std::setprecision(16)<<std::fixed;
			}

		}

	cout<<"Num Forest edges: "<<num_edges(graphFOF)<<endl;
	std::vector<int> component(num_vertices(graphFOF));
	int num = connected_components(graphFOF, &component[0]);
	cout<<"Number of components:"<<num<<endl;
	cout<<"======================="<<endl;
	/////////////// BUILD The MST tree
	weight_map_type weight = get(edge_weight, graphFOF);
	std::vector < Edge > spanning_tree;
	Graph graphMST;
	std::cout << "Building MST:." ;	
	kruskal_minimum_spanning_tree(graphFOF, std::back_inserter(spanning_tree));
	std::cout <<"MST num edges="<< spanning_tree.size()<<" .. done " << std::endl;
	cout<<"Make graphMST from Spanning Edges"<<endl;
	for (std::vector < Edge >::iterator ei = spanning_tree.begin();
			ei != spanning_tree.end(); ++ei) {
				//if(get(edge_weight, graphFOF, *ei)<1.0)
			add_edge(source(*ei, graphFOF),target(*ei, graphFOF),get(edge_weight, graphFOF, *ei),graphMST);
//			fout<<source(*ei, graphFOF)<<" "<<target(*ei, graphFOF)<<" "<<get(edge_weight, graphFOF, *ei)<<endl;
		}
	fout.close();
	/////////////////////////////////////////////////////////
	cout<<"Getting connected components"<<endl;
	num = connected_components(graphMST, &component[0]);
	cout<<"After MST: number of components:"<<num<<endl;
	cout<<"======================="<<endl;
exit(0);

	}
///////////////////////////////////////////////////////////////////
void read_ascii_n(string fname, array2dfloat &data, uint dim, int &ip)
	{
	std::ifstream fin(fname.c_str());
	unsigned int i=0,np;
	float tf;
	ip=0;
	fin>>np;
	data.resize(boost::extents[np][dim]);
	for(uint il=0;il<np;il++)
		{	
		for(i=0;i<dim;i++)
			{fin>>tf;data[ip][i]=tf;}
		ip++;
		}
	fin.close();
	cout<<"Got Np= "<<data.size()<<endl;
	}
#include <boost/format.hpp>
void dumpNgbTofile(const std::string file, const int *pNGB,const int NUM_NGB)
	{	
	using boost::format;
	using boost::io::group;
	std::ofstream of(file.c_str());
	for(int i=0;i<NUM_NGB;i++)
		of<< pNGB[i]<< "  ";
	of.close();
	
	}
using namespace boost;
	typedef adjacency_list < vecS, vecS, undirectedS,
		no_property, property < edge_weight_t, MyFloat > > Graph;
	typedef property_map<Graph, edge_weight_t>::type EdgeWeightMap;
	typedef property_map<Graph, edge_weight_t>::type	weight_map_type;
	typedef graph_traits < Graph >::edge_descriptor Edge;
	typedef graph_traits < Graph >::vertex_descriptor Vertex;
	typedef std::pair<int, int> E;
	
	void DumpEdges(Graph &g, std::string name="c:/arm2arm/DATA/edge_test.txt")
	{
	std::ofstream fout(name.c_str());
	//run over all edges and dump
	graph_traits<Graph>::edge_iterator eiter, eiter_end;
	for (tie(eiter, eiter_end) = edges(g); eiter != eiter_end; ++eiter) {
		fout << source(*eiter, g) << "  " << target(*eiter, g)<<" "<<get(edge_weight, g, *eiter)<<endl;
		}
	fout.close();
	};

template <typename EdgeWeightMap>
struct less_edge_weight {
  less_edge_weight() { }
  less_edge_weight(EdgeWeightMap weight,MyFloat thr=5.0) : m_weight(weight),m_threshold(thr) { }
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return m_threshold > get(m_weight, e);
  }
  EdgeWeightMap m_weight;
  MyFloat m_threshold;
};


void init_test_data()
	{
	cout<<"Starting to build FOF Graph based on adaptive Ngb list"<<endl;
	Graph graphFOF;
	int ntest=10;
	int i;
	array2dfloat realdata;
	read_ascii_n("c:/arm2arm/DATA/test_mahal.txt",realdata,ND_GROUPS,ntest); 
	
		All.DesNumNgb=3;
		All.NumPart=ntest;
		All.MaxPart=ntest;
		allocate_memory();
		for(i=0;i<ntest;i++)
			{
			 Part[i].id=i;
			 memset(&Part[i].Pos[0], 0, sizeof(Part[i].Pos[0])*6);
			}
		if(true)
			{
			for(int i=0;i<ntest;i++)
				{
				for(int j=0;j<ND_GROUPS; j++)
					{
					 Part[i].Pos[j]=(float)realdata[i][j];
					}
				}
			}else{
		Part[0].Pos[0]=-1.3557630f;Part[0].Pos[1]=-47.625917f;
		Part[1].Pos[0]=-1.7062211f;Part[1].Pos[1]=-63.574198f;
		Part[2].Pos[0]=-1.8655050f;Part[2].Pos[1]=-75.086588f;
		Part[3].Pos[0]=-1.0181101f;Part[3].Pos[1]=-16.815058f;
		Part[4].Pos[0]=1.3764034f;Part[4].Pos[1]=53.766238f;
		Part[5].Pos[0]=2.4233045f;Part[5].Pos[1]=49.306687f;
		Part[6].Pos[0]=2.1975691f;Part[6].Pos[1]=61.207182f;
		Part[7].Pos[0]=0.99472014f;Part[7].Pos[1]=60.852764f;
		Part[8].Pos[0]=-2.1623941f;Part[8].Pos[1]=-46.242952f;
		Part[9].Pos[0]=1.1159960f;Part[9].Pos[1]=24.211842f;
			}
		int dim=ND_GROUPS;
		realdata.resize(boost::extents[All.NumPart][dim]);
		for (i=0; i<All.NumPart; i++) {
			for (int j=0; j<dim; j++) 
				realdata[i][j] = Part[i].Pos[j];
			}
		std::vector<MyFloat> qv;
		int iq;
		const int NUM_NGB=15;
		kdtree2_result_vector ngblist, ngblistbrut;
		std::cout << "Build the KD-Tree:.." ;
		kdtree2*  tree = new kdtree2(realdata);
		tree->sort_results=true;

		boost::numeric::ublas::matrix<float> ngbGroup(ND_GROUPS,NUM_NGB);	
		CKernel<float> *pKernel=new CEpanechikov<float>(ND_GROUPS);
		int j=0;
		for(i=0;i<All.NumPart;i++)
			{
			iq=i;
			qv.clear();
			qv.push_back(Part[iq].Pos[0]);
			qv.push_back(Part[iq].Pos[1]);
			Part[i].pNGB=new int[NUM_NGB];
			tree->n_nearest( qv , NUM_NGB,ngblist);		
			int ingb=0;
			BOOST_FOREACH( kdtree2_result ngbNum, ngblist )
				{	
				j=ngbNum.idx;
				Part[i].pNGB[ingb++]=j;
				}
			//fill_ngb_group<float>(i,ngbGroup);
			//Metrics::MahalDistance<float> Mahal(ngbGroup);
			
			float  dist=0.0;
			
			//cout<<i<<endl;
			//cout<<"\t\t";
			float hsml2=(float)ngblist[NUM_NGB-1].dis;//Mahal.doDist(&Part[i].Pos[0],&Part[ngblist[NUM_NGB-1].idx].Pos[0]);
			float hsml=std::sqrt(hsml2);
			float u=0.0f, mahalDet=0.0f;
			Part[i].Est=0.0;
			BOOST_FOREACH( kdtree2_result ngbNum, ngblist )
				{	
				j=ngbNum.idx;
				//cout<<ngbNum.idx<<" ("<<ngbNum.dis<<") ";
				dist=ngbNum.dis;//Mahal.doDist(&Part[i].Pos[0],&Part[j].Pos[0]);
				Part[i].Est+=(pKernel->W(std::sqrt(dist/hsml2)));
				add_edge(i,ngbNum.idx,dist,graphFOF);
				mahalDet+=ngbNum.metric->getDet();
				}
			float hd=1.0;
			for(int idim=0;idim<ND_GROUPS;idim++)
				{
				hd*=hsml;
				}
				Part[i].Est/=(hd*std::sqrt(mahalDet/(float)NUM_NGB));
			
				if(i==5)dumpNgbTofile("c:/arm2arm/DATA/ngb.txt", &Part[i].pNGB[0], NUM_NGB);
			//cout<<endl;
			}

/////////////////////////////////////////
	cout<<"Num Forest edges: "<<num_edges(graphFOF)<<endl;
	std::vector<int> component(num_vertices(graphFOF));
	int num = connected_components(graphFOF, &component[0]);
	cout<<"Number of components:"<<num<<endl;
	cout<<"======================="<<endl;
	DumpEdges(graphFOF);
	/////////////////////////////////
/////////////// BUILD The MST tree
	weight_map_type weight = get(edge_weight, graphFOF);
	std::vector < Edge > spanning_tree;
	Graph graphMST;
	std::cout << "Building MST:." ;	
	kruskal_minimum_spanning_tree(graphFOF, 
								  std::back_inserter(spanning_tree)
								  );
	std::cout <<"MST num edges="<< spanning_tree.size()<<" .. done " << std::endl;
	cout<<"Make graphMST from Spanning Edges"<<endl;
	
	
	for (std::vector < Edge >::iterator ei = spanning_tree.begin();
			ei != spanning_tree.end(); ++ei) {
				//if(get(edge_weight, graphFOF, *ei)<1.0)
			add_edge(source(*ei, graphFOF),target(*ei, graphFOF),get(edge_weight, graphFOF, *ei),graphMST);
//			fout<<source(*ei, graphFOF)<<" "<<target(*ei, graphFOF)<<" "<<get(edge_weight, graphFOF, *ei)<<endl;
		}
	DumpEdges(graphMST, "c:/arm2arm/DATA/edge_test_mst.txt");
//////////////// NEW version for connected components
	less_edge_weight<EdgeWeightMap> filter(get(edge_weight, graphMST), 200.0);
	filtered_graph<Graph, less_edge_weight<EdgeWeightMap> > fg(graphMST, filter);
	component.resize(num_vertices(fg));
	cout<<"Filtered grapth has a: "<<connected_components(fg, &component[0])<<endl;
	}
///////////////////////////////////////////////////////////////////
//"snap_gal_sfr_0450.ascii"	
int main(int argc,char **argv) {

	// Read and check command line parameters.
	cimg_usage("Compute a HOP over the particles with given Est and Rho files");
	const char *file_ascii = cimg_option("-ascii_file",(char*)0,"Input in Ascii format");
	const char *file_gadget2;
	if( SMALL_TEST)
		file_gadget2= cimg_option("-gadget2_file","C:\\arm2arm\\DATA\\MODEL7\\MODELS\\MODEL7\\RUNG2\\SNAPS\\test_0450","Input in Gadget-2 format");
	else
		file_gadget2= cimg_option("-gadget2_file","C:\\arm2arm\\DATA\\MODEL7\\MODELS\\MODEL7\\RUNG2\\SNAPS\\snap_gal_sfr_0450","Input in Gadget-2 format");
	const char *file_base  = cimg_option("-base4files","C:/arm2arm/DATA/MODEL7/MST_GRAPH/","Base for input and output");
	const char *HOP_file  = cimg_option("-hopfile","hop_0450","HOP algorithm output,eg: groups and their IDs.");
	const bool visu     = cimg_option("-visu",true,"Visualization mode");
	const int shape  = cimg_option("-s",1,"shape [0,6]");
	const int profile  = cimg_option("-p",0,"profile [0,7]");

	////////////////////////////////////////////////////////////////////
	typeVecData data;

	base=string(file_base);
	string hopfile=base+string("/")+string(HOP_file);
	
	const int dim=3;
	
	init_test_data();
exit(0);
	if(file_ascii==0)
		read_gadget(file_gadget2, &data);

	////////////////////////////////////////////////
	// FOF_6D by Mahalanobis distance;
	
	MahalanobisFOF();
	exit(0);
	////////////////////////////////////////////////
	/// Smooth Est with 64 Ngb
	SmoothSph();
	/////////////////////// GET setAB This is the HOP stuff based on Enbid Density  one can test also for RHO by SPH//////////////////////////////////////

	MyFloat alpha=0.05;
	CHOP *hop=new CHOP(alpha, MIN_NGRP);

	int seeds=hop->get_seeds();
	cout<<"We got "<<seeds<<" seeds out of "<< All.NumPart<<endl;


	hop->write_catalogues(hopfile);

	///////////////////////////////////////////////////////////////
	// Here we need to cut and store the catalogues
	///////////////////////////////////////////////////////


	////////////////////////////////
	return EXIT_SUCCESS;
	}