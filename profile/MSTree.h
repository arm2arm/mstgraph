#ifndef _MY_MSTREE_
#define _MY_MSTREE_
#include "kdtree2.hpp"
#include <iostream>
#include <string>
#include <vector>
#include "utils.h"
using std::cout;
using std::endl;
using std::cerr;
using std::vector;

#include <boost/timer.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
///////////////////
#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/edge_connectivity.hpp>
#include <boost/graph/filtered_graph.hpp>
#include "MSTGroup.h"
///////////////////
#define MyFloat double
using namespace boost;
typedef adjacency_list < vecS, vecS, undirectedS,
no_property, property < edge_weight_t, MyFloat > > Graph;
typedef property_map<Graph, edge_weight_t>::type EdgeWeightMap;
typedef property_map<Graph, edge_weight_t>::type	weight_map_type;


typedef graph_traits < Graph >::edge_descriptor Edge;
typedef graph_traits < Graph >::vertex_descriptor Vertex;
typedef std::pair<int, int> E;

///////////////////

class CMSTree
	{
	typedef boost::multi_array<MyFloat,2> array2dfloat;
	public:
		//CMSTree(void);
		CMSTree(vector<float> &x, vector<float> &y,vector<float>  &z,
			float afof_eps=0.05,
			size_t min_num=100,
			size_t MAXNGB=15, bool update=true, std::string fname_="mst_group"):
		m_x(x),m_y(y),m_z(z),dim(3),m_afof_eps(afof_eps),m_min_num(min_num),m_maxNGB(MAXNGB), 
			fname(fname_),
			tree(NULL), m_verbose(false){
				fnameidx=fname+".grp.idx";
				fname+=".grp.txt";
				//if(!update && !load())
					{
					cout<<"geting NGB:"<<m_maxNGB<<endl;
					return;
					FillData();
					BuildKDTree();
					BuildGraph();
					GetMST();
					//CompileCATS();
					}
			};
		~CMSTree(void);
		void GetMST()
			{
			scoped_timer timemme("Get MSTGraph:.....");
			/////////////// BUILD The MST tree
			std::vector < Edge > spanning_tree;
			std::cout << "Building MST:." ;	
			kruskal_minimum_spanning_tree(graphFOF, std::back_inserter(spanning_tree));
			std::cout <<"MST num edges="<< spanning_tree.size()<<" .. done " << std::endl;

			int N = spanning_tree.size();
			typedef adjacency_list<vecS, vecS, undirectedS> UndirectedGraph;
			UndirectedGraph g(N);
			for(int i=0;i<N;i++)
				add_edge(spanning_tree[i].m_source,spanning_tree[i].m_target, g);

			std::vector<int> component(num_vertices(g));
			size_t num = connected_components(g, &component[0]);

			std::vector<int>::size_type i;
			cout << "Total number of components: " << num << endl;
			m_MSTCatalog.resize(num);
			for ( i = 0; i != component.size(); ++i)
				{
				//cout << "Vertex " << i <<" is in component " << component[i] << endl;
				m_MSTCatalog[component[i]].insert(i,m_x[i],m_y[i],m_z[i], rho[i]*rho[i] );				
				}

			cout<<"done compiling"<<endl;


			//sort them by size
			std::sort(m_MSTCatalog.begin(),m_MSTCatalog.end()); 
			num = std::count_if(m_MSTCatalog.begin(), m_MSTCatalog.end(), IfGt<int>(m_min_num));
			cout << "Total number of components where Np>"<<m_min_num<<" : " << num << endl;
			for (i = 0; i < num; ++i)
				{// compile position and velocity
				m_MSTCatalog[i].DoneInsert(i);
				if(m_verbose)
					cout<<m_MSTCatalog[i]<<endl;
				}
			if(!m_verbose && m_MSTCatalog.size()!=0)
				cout<<m_MSTCatalog[0]<<endl;
			cout <<"Done catalogues"<< endl;


			};
		void FillData()
			{
			size_t i;
			N=m_x.size();
			realdata.resize(boost::extents[N][dim]);
			for (i=0; i<N; i++) {				
				realdata[i][0] = m_x[i];
				realdata[i][1] = m_y[i];
				realdata[i][2] = m_z[i];				
				}

			}
		void BuildKDTree()
			{
			std::cout << "Build the KD-Tree:.." ;
			tree = new kdtree2(realdata);
			tree->sort_results = true;
			std::cout<<".DONE."<< std::endl;
			}
		template<class TF>
		TF Wsph(TF rr, TF h)
			{
			TF retval=0;
			TF u=rr/h;
			if(0<=u&&u<=1 )
				retval=1.0f-3.0f/2.0f*u*u+3.0f/4.0f*u*u*u;
			else
				if(1<u&&u<=2)
					retval=1.0f/4.0f*(2.0f-u)*(2.0f-u)*(2.0f-u);
				else
					return 0;
			return retval/(3.1456f*h*h*h);
			}
		void BuildGraph()
			{
			std::cout << "Build the FOFGraph..." ;
				{
				scoped_timer timemme("Get FOFGraph:.....");
				std::vector<MyFloat> qv;
				size_t i, j;
				float  m_afof_eps2=m_afof_eps*m_afof_eps;
				vector<bool> is_done(N,false);
				std::cout << "Build Graph:" << realdata.size()<<" particles"<<"...";cout.flush();
				qv.resize(3);
				//	CKernel<float> *pKernel = new CEpanechikov<float>(3); // Init Kernel for calculations
				rho.resize(N,0);
				std::fill(rho.begin(), rho.end(), 0.0f);
				cout<<endl;
				size_t imax=std::min<size_t>((size_t)m_maxNGB, (size_t)5);
				for(i=0;i<N;i++){
					qv[0]=m_x[i];qv[1]=m_y[i];qv[2]=m_z[i];
					tree->n_nearest( qv, m_maxNGB, ngblist);
					for(size_t ingb=1;ingb<imax;ingb++)
						{
						j=ngblist[ingb].idx;
						if(ngblist[ingb].dis<=m_afof_eps2)
							add_edge(i,j,graphFOF);
						}
					double hsml=sqrt(ngblist[ngblist.size()-1].dis);
					for(size_t ingb=0;ingb<ngblist.size();ingb++)
						{
						rho[i]+=(float)
							Wsph<double>(sqrt(ngblist[ingb].dis) , hsml);

						}
					if(m_verbose)
						if(i%100==0)cout<<i<<"\r";
					}

				if(m_verbose)
					report_components( graphFOF,false);
				//delete pKernel;
				std::cout<<".DONE."<< std::endl;
				}
			}

		float Distance2(size_t i, size_t j)
			{
			float distx=(m_x[i]-m_x[j]);
			float disty=(m_y[i]-m_y[j]);
			float distz=(m_z[i]-m_z[j]);
			return distx*distx+disty*disty+distz*distz;
			}
		template <class Graph>
		void report_components(Graph& g, bool full_verbose=false) 
			{
			if(full_verbose)
				scoped_timer timemme("Get Components:.....");
			std::vector<int> component(num_vertices(g));
			int num = connected_components(g, &component[0]);

			std::vector<int>::size_type i;
			cout << "Total number of components: " << num << endl;
			if(full_verbose)
				{
				for (i = 0; i != component.size(); ++i)
					cout << "Vertex " << i <<" is in component " << component[i] << endl;
				}
			cout << endl;

			}
		void dump(int ig)
			{
			std::string fname=std::string("part_ig")+boost::lexical_cast<std::string>(ig)+std::string(".ascii");
			std::ofstream of(fname.c_str());
			if(of.is_open())
				for(size_t i=0;i<m_MSTCatalog[ig].id.size();i++)
					{
					int ip=m_MSTCatalog[ig].id[i];

						{
						of<<std::setw(32)<<std::setprecision(16)<<m_x[ip]<<" "<<m_y[ip]<<" "<<m_z[ip]<<" "<<rho[ip]<<endl;
						}
					}
				of.close();
			}


		///////////////////////
		void save(){

			std::ofstream of(fname.c_str());
			std::ofstream ofidx(fnameidx.c_str());

			if(of.is_open())
				{
				if(m_MSTCatalog.size()!=0)
					{
					of<<m_MSTCatalog[0]<<endl;
					std::copy(m_MSTCatalog[0].id.begin(),
						m_MSTCatalog[0].id.end(), std::ostream_iterator<double>(ofidx," "));
					}
				}

			of.close();
			ofidx.close();
			}
		bool load(){
			bool result=true;
			////
			std::ifstream stream(fname.c_str());
			std::ifstream streamidx(fnameidx.c_str());
			if(!stream.is_open()){
				cerr<<"\nCannot open group file: "<<fname<<endl;
				return false;
				}
			if(!stream.is_open()){
				cerr<<"\nCannot open indexes for group file: "<<fname<<endl;
				return false;
				}
			CMSTGroup grp; 
			unsigned int idx;
			stream>>grp;
			for(int i=0;i<grp.Ntotal;i++)
				{streamidx>>idx;grp.id.push_back(idx);}
			m_MSTCatalog.push_back(grp);
			stream.close();
			streamidx.close();
			////
			return result;
			}
		///////////////////////

	protected:
		/////////////////////////
		int m_maxNGB;
		kdtree2_result_vector ngblist, ngblistbrut;
		array2dfloat realdata; 
		size_t N;
		size_t dim;
		float m_afof_eps;
		size_t m_min_num;
		////////////////
		kdtree2*  tree;
		Graph graphFOF;
		bool m_verbose;
		std::string fname;				    
		std::string fnameidx;

	public:
		vector<float> &m_x;
		vector<float> &m_y;
		vector<float> &m_z;
		vector<float> rho;
		TMSTCat m_MSTCatalog;


	};


#endif
