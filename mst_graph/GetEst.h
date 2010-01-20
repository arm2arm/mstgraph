#ifndef GET_EST
#define GET_EST
#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>



#include <boost/numeric/ublas/matrix.hpp>
#include <cstring>
#include "ANN/ANN.h"
#include "functions.h"
#include "program_settings.h"
#include "readers.h"
//extern enum eDOSORT;

//enum eDOSORT{BY_RHO, BY_EST,BY_AEST, BY_POS};

using namespace boost;
typedef adjacency_list < boost::vecS, vecS, undirectedS,
no_property, property < edge_weight_t, float > > Graph;
typedef property_map<Graph, edge_weight_t>::type EdgeWeightMap;
typedef property_map<Graph, edge_weight_t>::type	weight_map_type;

typedef graph_traits < Graph >::edge_descriptor Edge;
typedef graph_traits < Graph >::vertex_descriptor Vertex;
typedef std::pair<int, int> E;

class GetEst
	{
	public:
		typedef struct{float Pos[3];} tARR;
		std::vector<tARR> ARR;
		GetEst();
		~GetEst(void);
		void run();
		void run_byKdTree();
		void run_graph_Ap();
		void run_ANN(std::vector<float> &annRho,int dim=6,//Dimensions
			int k=40,//Estimation Neighbours
			int bucketsize=4,//Tuning parameter for ANN speed
			bool verbose=false
			);

		void Run_SPHEst();
		void SmoothByANN(int dim,std::vector<float> &annEst,//What to smooth
						 std::vector<float> &annSmooth,//result
						 int nsph,//Smoothing Neighbours in 3D
						 bool verbose=true
						 );

		template<class T>
		void fill_ngb_group(int i,boost::numeric::ublas::matrix<T> &ngbGroup, int count , eDOSORT which_est=BY_EST);
		template<class T>
		void fill_ngb_group_Ap(int i,boost::numeric::ublas::matrix<T> &ngbGroup,int count, eDOSORT which_est);

		void DumpVectorEst(std::string filename, int which=0);
		template <class T>
		void DumpVector(std::string fname, std::vector<T> &vec);
		template <class T>
		void LoadDumpVector(std::string fname, std::vector<T> &vec);
		void run_graph();
		void makeMSTree();
		void DumpNgb(std::string fname);
		void LoadNgb(std::string file);
		void AssignOneParticle(std::vector< std::vector<double> > &image,CKernel<double> *pKernel, 
							   double X, double Y, double Z, double Rho,double Hsml);
		void PutOnGrid(std::vector< std::vector<double> > &Grid2D, //Out 
					   std::vector<float> &annEst,//In: What to put on grid
					   std::vector<double> &Region,// In: CENT[3], R;//size=4 or size=6, eg: XC,YC,ZC,Rx, Ry, Rz;
					   // or if Size=3; Xc, YC, R
					   int nsph=0,//Smoothing Neighbours in 3D
					   bool verbose=true
					   );
		void SaveImage(std::string fname,std::vector< std::vector<double> > &Grid2D);
		void DumpMstCat(Graph &graphFOF,std::vector < Edge > &spanning_tree);
		int dim;
		static const int NUM_NGB=10;
		size_t Klink; //Ngb numbers used for linking;
		std::vector<float> Rho, myEst;
		std::vector<float> annRho, annHsml;
		std::vector<float> annEst;
		std::vector<float> annSmooth;

		typedef struct tagTNgb{std::vector< std::pair<int, float> > p;} TNgb;
		std::vector<TNgb> Mngb;
		/*! ANN structures
		ANNpointArray		dataPts;				// data points
		ANNpoint			queryPt;				// query point
		ANNidxArray			nnIdx;					// near neighbor indices
		ANNdistArray		dists;					// near neighbor distances
		ANNkd_tree*			kdTree;					// search structure

		*/
		ANNpointArray		dataPts;				// data points
		ANNpoint			queryPt;				// query point
		ANNidxArray			nnIdx;					// near neighbor indices
		ANNdistArray		dists;					// near neighbor distances
		ANNkd_tree*			kdTree;					// search structure
		void AllocateANNTreeStructures(int dim,int nPts, int kd);
	inline	void DeallocateANN(){	
		if(nnIdx!=NULL)
			{
			delete [] nnIdx;							// clean things up
			nnIdx=NULL;
			delete [] dists;
			delete kdTree;
			annClose();	
			}// done with ANN
			};
	void FillANNData(int dim);
	
	};



#endif
