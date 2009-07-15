//=======================================================================
// Copyright 2009 Arman Khalatyan, OAMP/LAM  
// distributed under GPL license
//=======================================================================


#include "mst_graph.h"
#include "ranker.h"
#include "coord.h"
#include "MSTGroup.h"
#include "kdtree2.hpp"


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
	bool operator()(ED w) const { return m_weights[w] > m_threshold; }

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

MyFloat CorrdDist(CCoord v1,CCoord v2)
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
		if(rand()/MyFloat(RAND_MAX) < 1.0)
		(*data).push_back( t);
		i++;
		}
	fin.close();
	fwin.close();
	cout<<"Got Np= "<<(*data).size()<<endl;
	}

/////////////////////////////
int rmain()
	{

	Graph graphFOF;
	const int num_nodes = 10000;
	typeVecData data;//(num_nodes);
	vector<int> idata;
	int i,j;
	int N;
	const int dim=3;
	//some random data
	//srand(0);
	//generate(data.begin(),data.end(),gencoord);
	//
	read_ascii("test.txt", &data);
	//read_ascii("testvel.txt", &data);
	//read_ascii("RVpVt.txt", &data);
	
	////////////////////////////////////////////////
	array2dfloat realdata; 
	N=data.size();
	realdata.resize(extents[data.size()][3]);
	for (i=0; i<N; i++) {
		for (int j=0; j<dim; j++) 
			realdata[i][j] = data[i].pos[j];
		}
	kdtree2_result_vector ngblist, ngblistbrut;
	std::cout << "Build the KD-Tree:" << std::endl;
	kdtree2*  tree = new kdtree2(realdata);
	tree->sort_results = true;
	//////////////////////////////////////////////
	if(N<20)print_data(data);
	std::vector<MyFloat> qv;
	int iq=0,mynum_edges=0;
	std::cout << "Build Graph:" << std::endl;
	FILE *fout=fopen("testvol.txt", "w");
	double fac=1.0/double(N);	
	for(i=0;i<N;i++)
		{
			iq=i;
		//if(N<20)
//			cout<<"# ";
		qv.clear();
		qv.push_back(data[iq].pos[0]);
		qv.push_back(data[iq].pos[1]);
		qv.push_back(data[iq].pos[2]);
		tree->n_nearest( qv , 3,ngblist);
		MyFloat vol=0, vol1=0;
		//tree->n_nearest_brute_force(qv, 2, ngblistbrut);
		//if(ngblist.size()>0)
			{	
			BOOST_FOREACH( kdtree2_result ngbNum, ngblist )
				{
				j=ngbNum.idx;
				vol+=ngbNum.vol;
				add_edge(i,j,CorrdDist(data[i],data[j]),graphFOF);
				}
			//tree->SetIsActive(iq,false);
			
			//if(N<20)
				//cout<<mynum_edges<<" ngbsize("<<ngblist.size()<<") "<<iq<<" ===> "<<j<<endl;			
			}
			fprintf(fout,"%g %g %g %g\n", data[i].pos[0], data[i].pos[1], data[i].pos[2],ngblist.size()/vol); 
		printf("%d \t %f5.3\r", i,i*fac);
		}
	fclose(fout);
///////////// NOW BUILD MST
	property_map < Graph, edge_weight_t >::type weight = get(edge_weight, graphFOF);
	std::vector < Edge > spanning_tree;
	std::cout << "Building MST:." ;	
	kruskal_minimum_spanning_tree(graphFOF, std::back_inserter(spanning_tree));
	std::cout << " .. done " << std::endl;
	/*for (std::vector < Edge >::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) {
			std::cout << source(*ei, graphFOF) << " <--> " << target(*ei, graphFOF)
				<< " with weight of " << weight[*ei]
			<< std::endl;
		}
		*/
	//////////////////////////////////////
	mynum_edges=spanning_tree.size();
	std::cout << "Graph:" << std::endl;
	std::cout <<" Num Edges: "<<mynum_edges<<std::endl;
	if(N<20)
		print_graph(graphFOF, get(vertex_index, graphFOF));
	std::cout << std::endl;
	if(N<20)
		dump_dot_file("test.dot", graphFOF);
	delete tree;
	
	if(N<20)
		myprint(graphFOF);
	print_graph_stats(graphFOF);
	/////////////////////////////////////////////////
	//	Cutting Graph
	//First Get Cuts len
	//run over all edges 
	int Nbin=10;
	vector<int> hist(Nbin, 0);
	vector<MyFloat> we;
	MyFloat lw=0.0;
	graph_traits<Graph>::edge_iterator eiter, eiter_end;
	for (tie(eiter, eiter_end) = edges(graphFOF); eiter != eiter_end; ++eiter) {
		  lw=std::log10(get(edge_weight, graphFOF, *eiter));
			we.push_back(lw);		  
			
		  }
	MyFloat mi=*min_element(we.begin(),we.end());
	MyFloat ma=*max_element(we.begin(),we.end());
	for(unsigned int i=0;i<we.size();i++)
		{
	     lw=we[i];
		   if(mi<lw)lw=mi;
		   if(ma>lw)lw=ma;
		   hist[(unsigned int)((lw-mi)/(ma-mi)*(Nbin-1))]++;	
		}
	//////////////////////////////////////////////////////////
	MyFloat meanW=(ma-mi)/0.5;
	{
	positive_edge_weight<EdgeWeightMap> filter(meanW, get(edge_weight, graphFOF));
		filtered_graph<Graph, positive_edge_weight<EdgeWeightMap> >
			fg(graphFOF, filter);

	std::vector<int> component(num_vertices(fg));
	int num = connected_components(fg, &component[0]);

	
	cout << "Total number of components: " << num <<" Wcut="<< meanW<< endl;

	Graph graphMST;
	for (std::vector < Edge >::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) {
			//std::cout << source(*ei, graphFOF) << " <--> " << target(*ei, graphFOF)
			//	<< " with weight of " << weight[*ei]
			//<< std::endl;
			add_edge(source(*ei, graphFOF),target(*ei, graphFOF),get(edge_weight, graphFOF, *ei),graphMST);
		}

	positive_edge_weight<EdgeWeightMap> filterMST(meanW, get(edge_weight, graphMST));
		filtered_graph<Graph, positive_edge_weight<EdgeWeightMap> >
			fgmst(graphMST, filterMST);

		
		string fname="test_"+toString(meanW)+".idx";
		cut_and_report(meanW,fgmst, fname);	
	dump_dot_file("msttest.dot", fgmst);
	}
 
	//////////////////////////////////////////////////////////
	 // Test real:
	//////////////////////////////////////////////////////////
	for(MyFloat dW=mi;dW<ma;dW+=(ma-mi)/10.0)
		{
		positive_edge_weight<EdgeWeightMap> filter(dW, get(edge_weight, graphFOF));
		filtered_graph<Graph, positive_edge_weight<EdgeWeightMap> >
			fg(graphFOF, filter);

		
		string fname="test_"+toString(dW)+".idx";
		cut_and_report(dW,fg, fname);
  	  
		if(N<20)  
		  dump_dot_file("test.dot", fg);

		}
/////////////////////////////////////////////////
//visual part
		{
	write_gts("test.gts", data, graphFOF);
	write_mst_gts("MSTtest.gts", data, graphFOF,spanning_tree);
		}
	return EXIT_SUCCESS;
	}




int main()
	{

	Graph graphFOF;
	typeVecData data;//(num_nodes);
	std::vector<MyFloat> fdata;
	int i,j;
	int N;
	const int dim=3;
	read_ascii_vecn("C:\\arm2arm\\DATA\\MODEL7\\MST_GRAPH\\test_450.ascii", 
		"C:\\arm2arm\\DATA\\MODEL7\\MST_GRAPH\\test_450.ascii_4_ph.est",&data);	
	
	////////////////////////////////////////////////
	array2dfloat realdata; 
	N=data.size();
	realdata.resize(extents[data.size()][3]);
	for (i=0; i<N; i++) {
		for (int j=0; j<dim; j++) 
			realdata[i][j] = data[i].pos[j];
		}
	kdtree2_result_vector ngblist, ngblistbrut;
	std::cout << "Build the KD-Tree:.." ;
	kdtree2*  tree = new kdtree2(realdata);
	tree->sort_results = true;
	cout<<".DONE."<< std::endl;
	//////////////////////////////////////////////
	std::vector<MyFloat> qv;
	int iq=0,mynum_edges=0;
	std::cout << "Build Graph:" << realdata.size()<<" particles"<<std::endl;
	for(i=0;i<N;i++)
		{
		iq=i;
		qv.clear();
		qv.push_back(data[iq].pos[0]);
		qv.push_back(data[iq].pos[1]);
		qv.push_back(data[iq].pos[2]);
		tree->n_nearest( qv , 10,ngblist);		
			{	
			BOOST_FOREACH( kdtree2_result ngbNum, ngblist )
				{
				j=ngbNum.idx;				
				add_edge(i,j,(data[i].w+data[j].w)*0.5,graphFOF);
				}	
			tree->SetIsActive(iq,false);
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
	
////////////////////////////////////
	//get quantile for Edge weigth
	const int Nbin=200;
	vector<MyFloat> WHist(Nbin+1), WHistX(Nbin+1);
	vector<MyFloat> tvec;
	
	for (std::vector < Edge >::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) {
			tvec.push_back(log10(weight[*ei]));
		}
	MyFloat ma=*max_element(tvec.begin(), tvec.end());
	MyFloat mi=*min_element(tvec.begin(), tvec.end());
	cout<<"WMin="<<mi<<"Wmax="<<ma<<endl;
	int c=tvec.size();
    for(uint i=0;i<tvec.size();i++)
		{
		 WHist[int((tvec[i]-mi)/(ma-mi)*Nbin)]+=1;		 		 		
		}
	for(uint i=0;i<Nbin;i++)
		WHistX[i]=i*(ma-mi)/MyFloat(Nbin)+mi;

		MyFloat Whistsum=accumulate(WHist.begin(), WHist.end(), 0.0);
		MyFloat psum=0.0;
		MyFloat aq[]={0.01,0.025,0.05,0.1,0.2,0.4};//PDF quantile in the fraction
		vector<MyFloat> WList;
		int aqsize=sizeof(aq)/sizeof(MyFloat);
		iq=0;
     	for(int j=WHist.size()-1;j>0;j--)
		{
		if((psum+=WHist[j]) >= Whistsum*aq[iq])
			{
			cout<<"iq="<<iq<<" j="<<j<<" "<<aq[iq]<<" "<<Whistsum*aq[iq]<<" "<<
				psum<<" "<<WHistX[j]<<" "<<WHist[j]<< endl;				
			WList.push_back(WHistX[j]);
			iq++;
			if(iq> aqsize-1)break;
			}
		}
	vector<MyFloat> taqvec( aq, &aq[ aqsize] );
    
   DumpVector("Whist.txt", WHist, WHist.size());
   DumpVector("WhistX.txt", WHistX, WHistX.size());
   DumpVector2("WList.txt", WList, taqvec,taqvec.size());
///////////////////////////////////////////////////////////////
// Here we need to cut and store the catalogues
	//////////////////////////////////////////////////////////
	 // Test real:
   //////////////////////////////////////////////////////////
// Brut force MST graph generation
   Graph graphMST;
   for (std::vector < Edge >::iterator ei = spanning_tree.begin();
	   ei != spanning_tree.end(); ++ei) {
		   add_edge(source(*ei, graphFOF),target(*ei, graphFOF),get(edge_weight, graphFOF, *ei),graphMST);
	   }
   component.resize(num_vertices(graphMST));
    num = connected_components(graphMST, &component[0]);
	cout<<"Number of components:"<<num<<endl;
cout<<"=========================="<<endl;
///////////////////////////////////////////////////////
typedef std::map<int , CMSTGroup> TMSTCat;
TMSTCat MSTCat;
///////////////////////////////////////////////////////
   for(uint i=0;i<WList.size();i++)
	   {
		   {
		   MyFloat meanW=pow(10.0,WList[i]);

		   Remover<weight_map_type>	r(get(edge_weight, graphMST), meanW);
		   remove_edge_if(r, graphMST);
	
		   uint numvert=num_vertices(graphMST);
		   std::vector<int> component(numvert);
		   int num = connected_components(graphMST, &component[0]);
			cout<<"==============\nTotal Vertexes="<<numvert<<endl;
		   cout << "Total number of connected components: " <<std::setprecision(5)<< num <<" Wcut="<< meanW<< endl;
//// make catalogues 	
		   for(int comp=0;comp<num;comp++)
			   {
			    MSTCat.insert(std::make_pair(comp, CMSTGroup(comp,&data)));
			   }
		   std::vector<int>::size_type i;
		   for (i = 0; i != component.size(); ++i)
			   {
			   MSTCat[component[i]].Insert(i);
			   }
			for(int comp=0;comp<num;comp++)
				{
				MSTCat[comp].DoneInsert();
				if(MSTCat[comp].Ntotal >10)
					cout<<MSTCat[comp];
				}

		   /*string fname="test_"+toString(meanW)+".idx";
		   cut_and_report(meanW,fgmst, fname);	
		   dump_dot_file("msttest.dot", fgmst);*/
///////////
		   }
	   }
   ////////////////////////////////
   return EXIT_SUCCESS;
	}
