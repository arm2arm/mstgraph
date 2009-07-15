#ifndef _MST_GRAPH_
#define _MST_GRAPH_


#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/timer.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>

#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/math/distributions.hpp>
//#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/error_of.hpp>
#include <boost/accumulators/statistics/error_of_mean.hpp>
//#include <boost/accumulators/statistics.hpp>
//#include <boost/accumulators/statistics/tail_quantile.hpp>

//using namespace boost::accumulators;
//typedef accumulator_set<double, stats<boost::accumulators::tag::tail_quantile<left> > > accumulator_t_left;
//typedef accumulator_set<double, stats<boost::accumulators::tag::tail_quantile<right> > > accumulator_t_right;

using namespace boost;
using namespace boost::math::policies;
using namespace boost::math;


#include <string>
#include <fstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <numeric>
#include <functional>
#include <iostream>
#include "cycle.h"


using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

typedef unsigned int uint;
typedef double MyFloat;

template<class T>
void DumpVector(std::string filename, std::vector<T> vec, unsigned int n=0,bool verbose=false)
	{
	std::ofstream of(filename.c_str());
	if(!of.good())assert("Error");
	unsigned int i=0, np=min(n,(unsigned int )vec.size());
	for(i=0;i<np;i++)
		{
		if(verbose)cout<<std::setprecision(2)<<std::fixed<<i<<") "<<vec[i]<<endl;
		of<<i<<" "<<vec[i]<<endl;
		}
	of.close();
	}
template<class T>
void DumpVector2(string filename, std::vector<T> vec,std::vector<T> vec1, unsigned int n=0,bool verbose=false)
	{
	std::ofstream of(filename.c_str());
	if(!of.good())assert("Error");
	unsigned int i=0, np=min(n,(unsigned int )vec.size());
	for(i=0;i<np;i++)
		{
		if(verbose)cout<<std::setprecision(2)<<std::fixed<<i<<") "<<vec[i]<<endl;
		of<<i<<" "<<vec[i]<<" "<<vec1[i]<<endl;
		}
	of.close();
	}

///////////////////////////////

class scoped_timer {
	boost::posix_time::ptime start_;
	std::string m_text;
public:    
	inline  void   SetText(std::string text){m_text=text;};
	scoped_timer(std::string text) 
		: start_(boost::posix_time::microsec_clock::universal_time()),m_text(text)
		{
		m_start= getticks();
		}
	~scoped_timer() {

		boost::posix_time::ptime stop( boost::posix_time::microsec_clock::universal_time() );
		m_end=getticks();
		double clock_cycles=elapsed(m_end, m_start);
		std::cout<<std::setprecision(2)<<std::fixed;
		std::cout <<" "<<m_text<< " done in " << ( stop - start_).total_milliseconds() << " milli seconds or "<<clock_cycles<<" CPU cycles "<<std::endl;
		}
protected:
	ticks m_start;
	ticks m_end;
	};


#endif