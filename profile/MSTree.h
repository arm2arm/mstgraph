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
#define MyFloat float
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
		CMSTree(vector<float> &x, vector<float> &y,vector<float>  &z,float afof_eps=0.05,size_t min_num=100,size_t MAXNGB=4):
		m_x(x),m_y(y),m_z(z),dim(3),m_afof_eps(afof_eps),m_min_num(min_num),m_maxNGB(MAXNGB), tree(NULL){			

			FillData();
			BuildKDTree();
			BuildGraph();
			GetMST();
			//CompileCATS();

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
			int num = connected_components(g, &component[0]);

			std::vector<int>::size_type i;
			cout << "Total number of components: " << num << endl;
			m_MSTCatalog.resize(num);
			for ( i = 0; i != component.size(); ++i)
				{
				//cout << "Vertex " << i <<" is in component " << component[i] << endl;
				m_MSTCatalog[component[i]].insert(i,m_x[i],m_y[i],m_z[i] );
				}
			// Lets remove small objects
			/*m_MSTCatalog.erase(
				std::remove_if(m_MSTCatalog.begin(),
				m_MSTCatalog.end(), IfLt<int>(m_min_num)
				)
				);*/
			cout<<"done compiling"<<endl;
			/* TMSTCat::iterator it = m_MSTCatalog.begin();
			while ( it != m_MSTCatalog.end() )
				{
				
				if ((*it).size() < m_min_num)
					{
					it = m_MSTCatalog.erase(it);
					}
				else
					{
					++it;
					}
				}
				*/
			cout << "Total number of components where Np>"<<m_min_num<<" : " << m_MSTCatalog.size() << endl;
//sort them by size
			std::sort(m_MSTCatalog.begin(),m_MSTCatalog.end()); 

			for (i = 0; i != m_MSTCatalog.size(); ++i)
				{// compile position and velocity
				m_MSTCatalog[i].DoneInsert();
				//cout<<m_MSTCatalog[i]<<endl;
				}
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
		void BuildGraph()
			{
			std::cout << "Build the FOFGraph..." ;
				{
				scoped_timer timemme("Get FOFGraph:.....");
				std::vector<MyFloat> qv;
				size_t i, j;
				float  m_afof_eps2=m_afof_eps*m_afof_eps;
				vector<bool> is_done(N,false);
				std::cout << "Build Graph:" << realdata.size()<<" particles"<<std::endl;
				qv.resize(3);
				for(i=0;i<N;i++){
					qv[0]=m_x[i];qv[1]=m_y[i];qv[2]=m_z[i];
					tree->n_nearest( qv, m_maxNGB, ngblist);
					for(size_t ingb=1;ingb<ngblist.size();ingb++)
						{
						j=ngblist[ingb].idx;
						if(ngblist[ingb].dis<=m_afof_eps2)
							add_edge(i,j,graphFOF);
						}
					}

				report_components( graphFOF,false);

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
		vector<float> &m_x;
		vector<float> &m_y;
		vector<float> &m_z;
		TMSTCat m_MSTCatalog;

	};


#endif