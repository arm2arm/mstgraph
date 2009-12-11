//=======================================================================
// Copyright 2009 Arman Khalatyan, OAMP/LAM  
// distributed under GPL license
//=======================================================================
#include <utility>
#include <map>
#include <set>
#include "mst_graph.h"
#include "ranker.h"
#include "coord.h"
#include "MSTGroup.h"
#include "kdtree2.hpp"
#include <math.h>
/////////////////////////////
#define MIN_NGRP 1000
#define MAXNGB 10
#define LL_LENGHT 0.1 // kpc maximum distance in real space
#define MAX_AP   70 // km/s maximum distance in the velocity space...
#define FAC_W   2.0

string base="C:\\arm2arm\\DATA\\MODEL7\\MST_GRAPH\\";
#define RAND_FRACTION 1.2
///////////////////////////////


typedef adjacency_list < vecS, vecS, undirectedS,
no_property, property < edge_weight_t, MyFloat > > Graph;
typedef property_map<Graph, edge_weight_t>::type EdgeWeightMap;
typedef property_map<Graph, edge_weight_t>::type	weight_map_type;


typedef graph_traits < Graph >::edge_descriptor Edge;
typedef graph_traits < Graph >::vertex_descriptor Vertex;
typedef std::pair<int, int> E;

///////////////////////////////
// SPH KERNEL///////
double Wsph(double rr, double h)
	{
	double retval=0;
	double u=rr/h;
	if(0<=u&&u<=1 )
		retval=1.0f-3.0f/2.0f*u*u+3.0f/4.0f*u*u*u;
	else
		if(1<u&&u<=2)
			retval=1.0f/4.0f*(2.0f-u)*(2.0f-u)*(2.0f-u);
		else
			return 0;
	return retval/(3.1456f*h*h*h);
	}
///////////////////////////////
///////////////////////////////
// For Cutter
template <typename EdgeWeightMap>
struct positive_edge_weight {
	positive_edge_weight() {}
	positive_edge_weight( MyFloat minW, EdgeWeightMap weight) : m_minW(minW),m_weight(weight) { }
	template <typename Edge>
	bool operator()(const Edge& e) const {
		MyFloat weight=get(m_weight, e);
		return m_minW > weight;
		}
	EdgeWeightMap m_weight;
	MyFloat m_minW;
	};


template<typename WMap>
class Remover
	{
	public:
		Remover(const WMap& weights, MyFloat threshold)
			: m_weights(weights), m_threshold(threshold) {}

		template<typename ED>
		bool operator()(ED w) const { return m_weights[w] < m_threshold; }

	private:
		const WMap&	m_weights;
		MyFloat			m_threshold;
	};

///////////////////////////////
template <typename T>
std::string toString(const T &thing) {
	std::ostringstream os;
	os << thing;
	return os.str();
	}





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
	/*CCoord loc_coord(counter,
	counter,
	counter,
	counter);*/
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
typedef multi_array<MyFloat,2> array2dfloat;

template <class Graph>
void myprint(Graph& g) {
	typename Graph::vertex_iterator i, end;
	typename Graph::out_edge_iterator ei, edge_end;
	for(boost::tie(i,end) = vertices(g); i != end; ++i) {
		cout << *i << " --> ";
		for (boost::tie(ei,edge_end) = out_edges(*i, g); ei != edge_end; ++ei)
			cout << target(*ei, g) << "  ";
		cout << endl;
		}
	}
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

/////////////////////////////
template <class Graph>
void cut_and_report(MyFloat dW, Graph& g, string fname, bool full_verbose=false) 
	{
	if(full_verbose)
		scoped_timer timemme("Get Components:.....");
	std::vector<int> component(num_vertices(g));
	int num = connected_components(g, &component[0]);

	std::vector<int>::size_type i;
	cout << "Total number of components: " << num <<" Wcut="<< dW<< endl;
	if(full_verbose)
		{
		for (i = 0; i != component.size(); ++i)
			cout << "Vertex " << i <<" is in component " << component[i] << endl;
		}
	cout << endl;

	}
///////////////////////////////
template <class Graph>
void print_graph_stats(Graph &g)
	{
	/*	typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
	for (edge_iterator e = edges(g).first; e != edges(g).second; ++e) {
	out << get(edge_weight, source(*e, g))<<"\n"; // HACK!
	}
	*/
	}
///////////////////////////////
template <class Graph>
void dump_dot_file(string fname, Graph &graphFOF)
	{
	if(true)
		{
		////////////////////////////////////////////////	
		// write
		dynamic_properties dp;
		dp.property("weight", get(edge_weight, graphFOF));
		dp.property("node_id", get(vertex_index, graphFOF));
		std::ofstream ofs( fname.c_str());
		std::cout << "write_graphviz:.....";
		write_graphviz(ofs, graphFOF, dp);
		ofs.close();
		std::cout << "Done" << std::endl;
		}else
			////////////////// WRITE BY HAND /////////////
		{
		///////////////////////////////////////////
		property_map < Graph, edge_weight_t >::type weight = get(edge_weight, graphFOF);
		std::vector < Edge > spanning_tree;

		kruskal_minimum_spanning_tree(graphFOF, std::back_inserter(spanning_tree));
		std::cout << "Print the edges in the MST:" << std::endl;
		for (std::vector < Edge >::iterator ei = spanning_tree.begin();
			ei != spanning_tree.end(); ++ei) {
				std::cout << source(*ei, graphFOF) << " <--> " << target(*ei, graphFOF)
					<< " with weight of " << weight[*ei]
				<< std::endl;
			}

		////////////////////////////////////////////////////////
		std::ofstream fout("kruskal.dot");
		fout << "graph A {\n"
			<< " rankdir=LR\n"
			<< " size=\"60,60\"\n"
			<< " ratio=\"filled\"\n"
			<< " edge[style=\"bold\"]\n" << " node[shape=\"circle\"]\n";
		graph_traits<Graph>::edge_iterator eiter, eiter_end;
		for (tie(eiter, eiter_end) = edges(graphFOF); eiter != eiter_end; ++eiter) {
			fout << source(*eiter, graphFOF) << " -- " << target(*eiter, graphFOF);
			if (std::find(spanning_tree.begin(), spanning_tree.end(), *eiter)
				!= spanning_tree.end())
				fout << "[color=\"blue\", label=\"" << get(edge_weight, graphFOF, *eiter)
				<< "\"];\n";
			else
				fout << "[color=\"red\", label=\"" << get(edge_weight, graphFOF, *eiter)
				<< "\"];\n";
			}
		fout << "}\n";
		fout.close();
			}
	}
////////////////////////////////////////
template <class Graph>
void write_gts(string fname, typeVecData &data, Graph &g)
	{
	std::ofstream fout(fname.c_str());
	fname+="_w";
	std::ofstream foutw(fname.c_str());
	unsigned int i,N=data.size();
	assert(N==num_vertices(g));
	fout<<num_vertices(g)<<" "<<num_edges(g)<<" "<<0<<
		" GtsSurface GtsFace GtsEdge GtsVertex "<<endl;
	for(i=0;i<N;i++)
		{
		fout<<data[i].pos[0]<<" "<<
			data[i].pos[1]<<" "<<
			data[i].pos[2]<<" "<<endl;
		}
	//run over all edges and dump
	graph_traits<Graph>::edge_iterator eiter, eiter_end;
	for (tie(eiter, eiter_end) = edges(g); eiter != eiter_end; ++eiter) {
		fout << source(*eiter, g) << "  " << target(*eiter, g)<<endl;
		foutw<< get(edge_weight, g, *eiter)<<endl;
		}

	fout.close();
	foutw.close();
	}

template <class Graph, class MST>
void write_mst_gts(string fname, typeVecData &data, Graph &g,MST &spanning_tree)
	{
	std::ofstream fout(fname.c_str());
	fname+="_w";
	std::ofstream foutw(fname.c_str());
	unsigned int i,N=data.size();
	//assert(N==num_vertices(g));
	fout<<N<<" "<<spanning_tree.size()<<" "<<0<<
		" GtsSurface GtsFace GtsEdge GtsVertex "<<endl;
	for(i=0;i<N;i++)
		{
		fout<<data[i].pos[0]<<" "<<
			data[i].pos[1]<<" "<<
			data[i].pos[2]<<" "<<endl;
		}
	//run over all edges and dump
	std::cout << "Print the edges in the MST:" << std::endl;
	property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
	for (MST::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) {			
			fout << source(*ei, g) << "  " << target(*ei, g)<<endl;
			foutw<< weight[*ei]<<endl;
		}

	fout.close();
	foutw.close();
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
bool compare_gt_MSTCat(const CMSTGroup &a, const CMSTGroup &b) 
	{
	return a.Ntotal > b.Ntotal;
	}

template <class T>
class NisGTThan{
public:
	NisGTThan(int incount):m_count(incount){};
	bool operator() (T a){
		return a.Ntotal>m_count;
		};
private:
	int m_count;

	};



///////////////////////////////////////////////////////////////////
void write_catalog(string catfile, TMSTCat *MSTCat)
	{

	TMSTCat::iterator::difference_type ResultOfCount=count_if((*MSTCat).begin(), (*MSTCat).end(),NisGTThan<CMSTGroup>(MIN_NGRP));
	cout<<"We Have TotalNgrp="<<(*MSTCat).size()<<endl;
	cout<<" But greather than "<<MIN_NGRP<<" = "<<ResultOfCount<<endl;
	std::ofstream fout(catfile.c_str());
	if(ResultOfCount<0)
		{
		cout<<"Not groups in the MST analysis"<<endl;
		return;
		}

	if(fout.bad()){cout<<"cannot open file for catalog:\n"<<catfile<<"\n"<<endl;exit(0);};
	uint i;
	for(i=0;i<(*MSTCat).size();i++)
		fout<<(*MSTCat)[i];

	fout.close();
	catfile+="_idx";
	fout.open(catfile.c_str(), std::ios::binary);

	if(fout.bad()){cout<<"cannot open file for catalog:\n"<<catfile<<"\n"<<endl;exit(0);};
	uint NGrp=ResultOfCount;//(*MSTCat).size();
	fout.write((char*)&NGrp,sizeof(NGrp));
	cout<<"We will write NGrp="<<NGrp<<" groups to file"<<endl;
	for(i=0;i<NGrp;i++)//(*MSTCat).size();i++)
		{
		fout.write((char*)&(*MSTCat)[i].Ntotal,sizeof(int));

		for(int j=0;j<(*MSTCat)[i].Ntotal;j++)
			fout.write((char*)&(*MSTCat)[i].id[j],sizeof(int));
		if(!(i%100))cout<<"iCat="<<i<<" Ntotal="<<(*MSTCat)[i].Ntotal<<"\r";
		}
	fout.close();
	}
////////////////////////////////////////////////////////////////////
template<class T> struct index_cmp {
index_cmp( T arr) : arr(arr) {}
bool operator()( size_t a,  size_t b) const
{ return arr[a] >arr[b]; }
 T arr;
};
///////////////////////////////////////////////////////////////////
typedef struct tagASTR{
	std::vector<std::pair<int,MyFloat>> setA;
	std::vector<int> setB;
	} TMyNgb;
vector<TMyNgb> NPart;

unsigned int count_seeds()
	{
	unsigned int i, counter=0;
	vector<uint> seed_vector;
	for(i=0;i<NPart.size();i++)
		{
			if(NPart[i].setB.size()==0)
				{
				seed_vector.push_back(i);
				counter++;

				}
		}
	DumpVector(base+string("/test_seeds.txt"), seed_vector, seed_vector.size());
	return counter;
	}
class TStr{
public:
	std::map<uint,uint> str;
	TStr(uint id,uint i){str.insert(std::pair<uint, uint>(id,i));};
	};

vector<unsigned int> sindex;
typedef  std::vector<MyFloat> TVMyFloat;
	TVMyFloat RhoVecNGB;

vector<TStr> cats;
void populate_structures(uint nstr)
	{
		
		uint i, id=1;
		MyFloat meanRho=(*std::max_element(RhoVecNGB.begin(),RhoVecNGB.end()))-(*std::min_element(RhoVecNGB.begin(),RhoVecNGB.end()));
		
		for(i=0;i<NPart.size();i++)
			{
			uint idx=sindex[i];// index sorted by Rho  
			if(NPart[idx].setB.size()==0)
				{
				TStr onestr(id,idx);
				cats.push_back(onestr);
				id++;
				}
			}
	
	}
	
int main()
	{

	Graph graphFOF;
	typeVecData data;//(num_nodes);
	std::vector<MyFloat> fdata;
	string snapfile=string("snap_gal_sfr_0450.ascii");
	int i,j;
	int N;
	MyFloat meanW=0.0;
	const int dim=3;
	string base="C:\\arm2arm\\DATA\\MODEL7\\MST_GRAPH\\";
	read_ascii_vecn(base+snapfile, 
		base+snapfile+string("_4_ph.est"),&data);	

	////////////////////////////////////////////////
	array2dfloat realdata; 
	N=data.size();
	MyFloat maxW=0.0;
	realdata.resize(extents[data.size()][3]);
	for (i=0; i<N; i++) {
		for (int j=0; j<dim; j++) 
			realdata[i][j] = data[i].pos[j];
		meanW+=data[i].w;
		maxW=max(data[i].w, maxW);
		}
	meanW/=(MyFloat)N;
	DumpVectorPos<CCoord>(base+string("/test_pos.txt"), data, data.size());
	kdtree2_result_vector ngblist, ngblistbrut;
	std::cout << "Build the KD-Tree:.." ;
	kdtree2*  tree = new kdtree2(realdata);
	tree->sort_results = true;
	cout<<".DONE."<< std::endl;
	//////////////////////////////////////////////
	std::vector<MyFloat> qv;
	int iq=0,mynum_edges=0;
	MyFloat dist=0.0;
	MyFloat Ap[4], dAp=0, dAr=0,  dW=0, dA=0.0;
	std::cout << "Build Graph:" << realdata.size()<<" particles"<<std::endl;
	for(i=0;i<N;i++)
		{

		iq=i;
		qv.clear();
		qv.push_back(data[iq].pos[0]);
		qv.push_back(data[iq].pos[1]);
		qv.push_back(data[iq].pos[2]);
		if(data[iq].w > meanW)
			{
			tree->n_nearest( qv , MAXNGB,ngblist);		
				{	
				BOOST_FOREACH( kdtree2_result ngbNum, ngblist )
					{
					j=ngbNum.idx;
					//dist=CoordDist(data[i],data[j]);
					//GetAP(data[i], &Ap[0]);
					//GetAP(data[j], &Ap[2]);

					//dAp=(Ap[1]-Ap[3]);
					//dAr=(Ap[0]-Ap[2]);
					//dA=sqrt(dAp*dAp+dAr*dAr);
					//dW=(data[i].w-data[j].w);
					//if( dA<=MAX_AP )// && dist<=LL_LENGHT)
					if(dist<=0.5)
					  add_edge(i,j, data[i].w,graphFOF);
					//dV=sqrt((data[i].vel[0]-data[j].vel[0])*(data[i].vel[0]-data[j].vel[0])+
					//(data[i].vel[1]-data[j].vel[1])*(data[i].vel[1]-data[j].vel[1])+
					//(data[i].vel[2]-data[j].vel[2])*(data[i].vel[2]-data[j].vel[2]));
					/*dW=(data[i].w-data[j].w);
					if( dist<LL_LENGHT && dW>0.00008)
					add_edge(i,j, dW,graphFOF);
					*/
					}	
				tree->SetIsActive(iq,false);
				}
			}
		if(i%10 ==0)
			{
			std::cout<<std::setprecision(2)<<std::fixed;
			std::cout<<i/double(N)*100.0<<"%"<<"\r";
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

	std::cout << "Building MST:." ;	
	kruskal_minimum_spanning_tree(graphFOF, std::back_inserter(spanning_tree));
	std::cout <<"MST num edges="<< spanning_tree.size()<<" .. done " << std::endl;
	/////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////
	// Here we need to cut and store the catalogues
	//////////////////////////////////////////////////////////
	// Test real:
	//////////////////////////////////////////////////////////
	// Brut force MST graph generation
	cout<<"Building MST from NGB groups.";std::cout.flush();
	Graph graphMST;
	MyFloat FAC_FORIT=1.1;
	MyFloat RhoThr=meanW*FAC_FORIT;
	uint It=0;
	num=0;
	while(num <2 && meanW <maxW)
		{
		for (std::vector < Edge >::iterator ei = spanning_tree.begin();
			ei != spanning_tree.end(); ++ei) {
				RhoThr=meanW*FAC_FORIT*It;
				if(get(edge_weight, graphFOF, *ei)>RhoThr)
					add_edge(source(*ei, graphFOF),target(*ei, graphFOF),get(edge_weight, graphFOF, *ei),graphMST);
			}
		cout<<"Done.."<<endl;
		cout<<"Get connected components..";std::cout.flush();
		component.resize(num_vertices(graphMST));
		num = connected_components(graphMST, &component[0]);
		cout<<"Number of components:"<<num<<endl;
		cout<<"..Done..It="<<It<<endl;
		It++;
		}

	
	cout<<"=========================="<<endl;
	///////////////////////////////////////////////////////
	for(uint iq=0;iq<1;iq++)
		{
			{
			MyFloat meanW=0.0001;

			//Remover<weight_map_type>	r(get(edge_weight, graphMST),meanW );
			//remove_edge_if(r, graphMST);

			uint numvert=num_vertices(graphMST);
			std::vector<int> component(numvert);
			int num = connected_components(graphMST, &component[0]);
			cout<<"==============\nTotal Vertexes="<<numvert<<endl;
			cout << "Total number of connected components: " <<std::setprecision(5)<< num <<" Wcut="<< meanW<< endl;
			//// make catalogues 	
			TMSTCat MSTCat(num);
			for(int comp=0;comp<num;comp++)
				{
				MSTCat[comp]= CMSTGroup(comp,&data);
				}
			std::vector<int>::size_type i;
			for (i = 0; i != component.size(); ++i)
				{
				MSTCat[component[i]].Insert(i);
				}
			for(int comp=0;comp<num ;comp++)
				{
				MSTCat[comp].DoneInsert();

				}
			cout<<"DONE"<<endl;
			cout<<"Sort Catalogues.";std::cout.flush();

			sort(MSTCat.begin(), MSTCat.end(), compare_gt_MSTCat);
			cout<<"DONE"<<endl;
			cout<<"Write up results..";std::cout.flush();

			string catfile=base+"MSTCat_6DFOF_450_4_00"+boost::lexical_cast<std::string>(iq);
			write_catalog(catfile, &MSTCat);
			cout<<"DONE"<<endl;
			/*string fname="test_"+toString(meanW)+".idx";
			cut_and_report(meanW,fgmst, fname);	
			dump_dot_file("msttest.dot", fgmst);*/
			///////////
			}
		}
	////////////////////////////////
	return EXIT_SUCCESS;
	}
