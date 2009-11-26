#ifndef __KDTREE2_HPP
#define __KDTREE2_HPP

//
// (c) Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004)
//
// Licensed under the Academic Free License version 1.1 found in file LICENSE
// with additional provisions in that same file.



//
// Implement a kd tree for fast searching of points in a fixed data base
// in k-dimensional Euclidean space.
//
// 
#include "mst_graph.h"
#include "metric_factory.hpp"
#include <vector>
#include <algorithm>

#define BOOST_NO_STDC_NAMESPACE

#include <boost/multi_array.hpp>
#include <boost/array.hpp>

typedef boost::multi_array<MyFloat,2>             kdtree2_array;
typedef boost::const_multi_array_ref<MyFloat,2>   kdtree2_ro_array;  // read only ref


typedef struct {
	MyFloat lower, upper;
	} interval;


// let the compiler know that this is a names of classes. 
class kdtree2_node; 
class searchrecord;

//
// struct KDTREE2_RESULT
// class  KDTREE2_RESULT
//


struct kdtree2_result {
	// 
	// the search routines return a (wrapped) vector
	// of these. 
	//

public:
	kdtree2_result():is_active(true){}
	MyFloat dis;  // its square Euclidean distance
	int idx;    // which neighbor was found
	bool is_active;
	MyFloat vol;
	Metrics::MahalDistance<MyFloat> *metric;
	}; 

class kdtree2_result_vector : public std::vector<kdtree2_result> {
	// inherit a vector<kdtree2_result>
	// but, optionally maintain it in heap form as a priority
	// queue.
public:

	//  
	// add one new element to the list of results, and
	// keep it in heap order.  To keep it in ordinary, as inserted,
	// order, then simply use push_back() as inherited
	// via vector<> 

	void push_element_and_heapify(kdtree2_result&);
	MyFloat replace_maxpri_elt_return_new_maxpri(kdtree2_result&);

	MyFloat max_value(); 
	// return the distance which has the maximum value of all on list, 
	// assuming that ALL insertions were made by
	// push_element_and_heapify() 
	};


//
// class KDTREE2
//
// The main data structure, one for each k-d tree, pointing
// to a tree of an indeterminate number of "kdtree2_node"s.
//

class kdtree2 {
public: 
	const kdtree2_array& the_data;   

	//pointing the dtata state is it active(true) or not(false);

	// "the_data" is a reference to the underlying multi_array of the
	// data to be included in the tree.
	//
	// NOTE: this structure does *NOT* own the storage underlying this.
	// Hence, it would be a very bad idea to change the underlying data
	// during use of the search facilities of this tree.
	// Also, the user must deallocate the memory underlying it.


	const int N;   // number of data points
	int dim; //
	bool sort_results;  // USERS set to 'true'. 
	const bool rearrange; // are we rearranging? 

public:
	//
	// constructor, passing in a multi_array<MyFloat,2> , aka
	// kdtree2_array. 
	//
	// constructor, has optional 'dim_in' feature, to use only
	// first 'dim_in' components for definition of nearest neighbors.
	//

	kdtree2(kdtree2_array& data_in, bool rearrange_in = false,int dim_in=-1);

	// destructor
	~kdtree2();

public:
	//remove i-th particle from NGB searcher.
	void SetIsActive(int i, bool state_in=false){data_state[i] = state_in;};
public:
	// search routines

	void n_nearest_brute_force(std::vector<MyFloat>& qv, int nn, kdtree2_result_vector& result);
	// search for n nearest to a given query vector 'qv' usin
	// exhaustive slow search.  For debugging, usually.

	void n_nearest(std::vector<MyFloat>& qv, int nn, kdtree2_result_vector& result);
	// search for n nearest to a given query vector 'qv'.

	void n_nearest_around_point(int idxin, int correltime, int nn,
		kdtree2_result_vector& result);
	// search for 'nn' nearest to point [idxin] of the input data, excluding
	// neighbors within correltime 

	void r_nearest(std::vector<MyFloat>& qv, MyFloat r2,kdtree2_result_vector& result); 
	// search for all neighbors in ball of size (square Euclidean distance)
	// r2.   Return number of neighbors in 'result.size()', 

	void r_nearest_around_point(int idxin, int correltime, MyFloat r2,
		kdtree2_result_vector& result);
	// like 'r_nearest', but around existing point, with decorrelation
	// interval. 

	int r_count(std::vector<MyFloat>& qv, MyFloat r2);
	// count number of neighbors within square distance r2.
	int r_count_around_point(int idxin, int correltime, MyFloat r2);
	// like r_count, c

	void kdtree2::n_density(std::vector<MyFloat>& qv, int nn, kdtree2_result_vector& result, MyFloat &den);
	// Get density for the qv point, taking NGB smoothing;

	friend class kdtree2_node;
	friend class searchrecord;
private:
	// private data members

	kdtree2_node* root; // the root pointer

	const kdtree2_array* data;
	// pointing either to the_data or an internal
	// rearranged data as necessary

	std::vector<int> ind; 
	std::vector<bool>  data_state; 
	// the index for the tree leaves.  Data in a leaf with bounds [l,u] are
	// in  'the_data[ind[l],*] to the_data[ind[u],*]

	kdtree2_array rearranged_data;  
	// if rearrange is true then this is the rearranged data storage. 


	static const int bucketsize = 100;  // 12 global constant. 

private:
	void get_metric(int l, int u, kdtree2_node* node);
	void set_data(kdtree2_array& din); 
	void build_tree(); // builds the tree.  Used upon construction. 
	std::ofstream fdump;
	kdtree2_node* build_tree_for_range(int l, int u, kdtree2_node* parent);
	void select_on_coordinate(int c, int k, int l, int u); 
	int select_on_coordinate_value(int c, MyFloat alpha, int l, int u); 
	void spread_in_coordinate(int c, int l, int u, interval& interv);
	};


//
// class KDTREE2_NODE
// 
// a node in the tree.  Many are created per tree dynamically.. 
//


class kdtree2_node {
public:
	// constructor
	kdtree2_node(int dim);
	//, int cut_dim_in,
	// 	       MyFloat cut_val_in, MyFloat cut_val_left_in, 
	//	       MyFloat cut_val_right_in);
	// destructor
	~kdtree2_node();
private:
	// visible to self and kdtree2.
	friend class kdtree2;  // allow kdtree2 to access private 
	int cut_dim;                                 // dimension to cut; 
	MyFloat cut_val, cut_val_left, cut_val_right;  //cut value
	int l,u;  // extents in index array for searching

	std::vector<interval> box; // [min,max] of the box enclosing all points
	Metrics::MahalDistance<MyFloat> metric;
	kdtree2_node *left, *right;  // pointers to left and right nodes. 

	void search(searchrecord& sr); 
	// recursive innermost core routine for searching.. 

	bool box_in_search_range(searchrecord& sr);
	// return true if the bounding box for this node is within the
	// search range given by the searchvector and maximum ballsize in 'sr'. 

	//void check_query_in_bound(searchrecord& sr); // debugging only

	// for processing final buckets. 
	void process_terminal_node(searchrecord& sr);
	void process_terminal_node_fixedball(searchrecord& sr);


	};



#endif
