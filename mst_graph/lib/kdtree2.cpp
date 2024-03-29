// Weight added in to the nodes to track 
// the density and the volume for the unit less multidimensional analysis
// Modified by A.Khalatyan 2009, LAM/OAMP, Marseille France


/////////////// ORIGINAL KD-Tree FROM ////////////////////////////////////////
// (c) Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004)
//
// Licensed under the Academic Free License version 1.1 found in file LICENSE
// with additional provisions in that same file.
///////////////////////////////////////////////////////////////////////////////

#include "kdtree2.hpp"

#include <algorithm> 
#include <cmath> 
#include <iostream>
#include <vector>
// utility
using std::cout;
using std::endl;
using std::max;
using std::min;

inline MyFloat squared(const MyFloat x) {
	return(x*x);
	}
template<class T>
inline void kdswap(T& a, T&b) {
	T tmp;
	tmp = a;
	a = b;
	b = tmp; 
	}




//
//       KDTREE2_RESULT implementation
// 

inline bool operator<(const kdtree2_result& e1, const kdtree2_result& e2) {
	return (e1.dis < e2.dis);
	}

//
//       KDTREE2_RESULT_VECTOR implementation
// 
MyFloat kdtree2_result_vector::max_value() {
	return( (*begin()).dis ); // very first element
	}

void kdtree2_result_vector::push_element_and_heapify(kdtree2_result& e) {
	push_back(e); // what a vector does.
	push_heap( begin(), end() ); // and now heapify it, with the new elt.
	}
MyFloat kdtree2_result_vector::replace_maxpri_elt_return_new_maxpri(kdtree2_result& e) {
	// remove the maximum priority element on the queue and replace it
	// with 'e', and return its priority.
	//
	// here, it means replacing the first element [0] with e, and re heapifying.

	pop_heap( begin(), end() );
	pop_back();
	push_back(e); // insert new
	push_heap(begin(), end() );  // and heapify. 
	return( (*this)[0].dis );
	}

//
//        KDTREE2 implementation
//

// constructor
kdtree2::kdtree2(kdtree2_array& data_in,bool rearrange_in,int dim_in)
: the_data(data_in),
N  ( data_in.shape()[0] ),
dim( data_in.shape()[1] ),
sort_results(false),
rearrange(rearrange_in), 
root(NULL),
data(NULL),
ind(N),
data_state(N)
	{
	//
	// initialize the constant references using this unusual C++
	// feature.
	//
	if (dim_in > 0) 
		dim = dim_in;

	build_tree();

	if (rearrange) {
		// if we have a rearranged tree.
		// allocate the memory for it. 
		printf("rearranging\n"); 
		rearranged_data.resize( boost::extents[N][dim] );

		// permute the data for it.
		for (int i=0; i<N; i++) {
			for (int j=0; j<dim; j++) {
				rearranged_data[i][j] = the_data[ind[i]][j];
				// wouldn't F90 be nice here? 
				}
			}
		data = &rearranged_data;
		} else {
			data = &the_data;
		}
	}


// destructor
kdtree2::~kdtree2() {
	delete root;
	}

// building routines
void kdtree2::build_tree() {
	for (int i=0; i<N; i++){ ind[i] = i;data_state[i]=true;} 
	fdump.open("c:/arm2arm/DATA/kd_dump.txt");
	root = build_tree_for_range(0,N-1,NULL);
	fdump.close();
	}

kdtree2_node* kdtree2::build_tree_for_range(int l, int u, kdtree2_node* parent) {
	// recursive function to build 
	kdtree2_node* node = new kdtree2_node(dim);
	// the newly created node. 

	if (u<l) {
		return(NULL); // no data in this node. 
		}


	if ((u-l) <= bucketsize) {
		// create a terminal node. 

		// always compute true bounding box for terminal node. 
		for (int i=0;i<dim;i++) {
			spread_in_coordinate(i,l,u,node->box[i]);
			}

		node->cut_dim = 0; 
		node->cut_val = 0.0;
		node->l = l;
		node->u = u;
		node->left = node->right = NULL;
		
		if(parent->childlist.size() >0)
			{
			node->metric=parent->metric;
			
			}else
				get_metric( l, u, node);
		
		} else {
			//
			// Compute an APPROXIMATE bounding box for this node.
			// if parent == NULL, then this is the root node, and 
			// we compute for all dimensions.
			// Otherwise, we copy the bounding box from the parent for
			// all coordinates except for the parent's cut dimension.  
			// That, we recompute ourself.
			//
			int c = -1;
			MyFloat maxspread = 0.0f;
			int m; 

			for (int i=0;i<dim;i++) {
				if ((parent == NULL) || (parent->cut_dim == i)) {
					spread_in_coordinate(i,l,u,node->box[i]);
					} else {
						node->box[i] = parent->box[i];
					}
				MyFloat spread = node->box[i].upper - node->box[i].lower; 
				if (spread>maxspread) {
					maxspread = spread;
					c=i; 
					}
				}

			// 
			// now, c is the identity of which coordinate has the greatest spread
			//

			if (false) {
				m = (l+u)/2;
				select_on_coordinate(c,m,l,u);  
				} else {
					MyFloat sum; 
					MyFloat average;

					if (true) {
						sum = 0.0f;
						for (int k=l; k <= u; k++) {
							sum += the_data[ind[k]][c];
							}
						average = sum / static_cast<MyFloat> (u-l+1);
						} else {
							// average of top and bottom nodes.
							average = (node->box[c].upper + node->box[c].lower)*0.5f; 
						}

					m = select_on_coordinate_value(c,average,l,u);
				}


			// move the indices around to cut on dim 'c'.
			node->cut_dim=c;
			node->l = l;
			node->u = u;
			//get_metric( l, u, node);// we dont calculate metric for all nodes.. only for terminal one generated once
			node->left = build_tree_for_range(l,m,node);
			node->right = build_tree_for_range(m+1,u,node);

			if (node->right == NULL) {
				for (int i=0; i<dim; i++) 
					node->box[i] = node->left->box[i]; 
				node->cut_val = node->left->box[c].upper;
				node->cut_val_left = node->cut_val_right = node->cut_val;
				} else if (node->left == NULL) {
					for (int i=0; i<dim; i++) 
						node->box[i] = node->right->box[i]; 
					node->cut_val =  node->right->box[c].upper;
					node->cut_val_left = node->cut_val_right = node->cut_val;
				} else {
					node->cut_val_right = node->right->box[c].lower;
					node->cut_val_left  = node->left->box[c].upper;
					node->cut_val = (node->cut_val_left + node->cut_val_right) / 2.0f; 
					//
					// now recompute true bounding box as union of subtree boxes.
					// This is now faster having built the tree, being logarithmic in
					// N, not linear as would be from naive method.
					//
					for (int i=0; i<dim; i++) {
						node->box[i].upper = max(node->left->box[i].upper,
							node->right->box[i].upper);

						node->box[i].lower = min(node->left->box[i].lower,
							node->right->box[i].lower);
						/*cout<<i<<" "<< node->box[i].upper<<" "<<
						node->box[i].lower<<" ";*/
						}
					}
		}

	//if( std::abs(((u-l)- bucketsize))<2 )
			{
			if(false)

				{fdump<<node->cut_dim<<" ";
			for (int i=0; i<dim; i++) {
				fdump<<"  "<<i<<" "<< node->box[i].lower<<" "<<
					node->box[i].upper<<" ";
				}
			if(node->left==NULL && node->right==NULL)fdump<<0;else fdump<<1;
			fdump<<endl;
				}
			}
	return(node);
	}

// this whould calculate adaptive metrics for a given node.
void kdtree2:: get_metric( int l, int u, kdtree2_node* node)
	{
	static uint times=0;
	int np=u-l;
	boost::numeric::ublas::matrix<MyFloat> ngbGroup(dim,np);
	//cout<<times<<")METRIC: "<<l<<" "<<u<<" on Np=" << np<<endl;
	for (int i=l,ip=0; i< u; i+=1) 
		{			 
		for(int k=0;k<dim;k++)
			ngbGroup(k,ip)=the_data[ind[i]] [k];
		node->childlist.push_back(ind[i]);
		ip++;
		}
	//cout<<ngbGroup<<endl;
	node->metric.Init(ngbGroup);
	times++;

	}
void kdtree2:: spread_in_coordinate(int c, int l, int u, interval& interv)
	{
	// return the minimum and maximum of the indexed data between l and u in
	// smin_out and smax_out.
	MyFloat smin, smax;
	MyFloat lmin, lmax;
	int i; 

	smin = the_data[ind[l]][c];
	smax = smin;


	// process two at a time.
	for (i=l+2; i<= u; i+=2) {
		lmin = the_data[ind[i-1]] [c];
		lmax = the_data[ind[i]  ] [c];

		if (lmin > lmax) {
			std::swap(lmin,lmax); 
			//      MyFloat t = lmin;
			//      lmin = lmax;
			//      lmax = t;
			}

		if (smin > lmin) smin = lmin;
		if (smax <lmax) smax = lmax;
		}
	// is there one more element? 
	if (i == u+1) {
		MyFloat last = the_data[ind[u]] [c];
		if (smin>last) smin = last;
		if (smax<last) smax = last;
		}
	interv.lower = smin;
	interv.upper = smax;
	//printf("Spread in coordinate %d=[%f,%f]\n",c,smin,smax);
	}


void kdtree2::select_on_coordinate(int c, int k, int l, int u) {
	//
	//  Move indices in ind[l..u] so that the elements in [l .. k] 
	//  are less than the [k+1..u] elmeents, viewed across dimension 'c'. 
	// 
	while (l < u) {
		int t = ind[l];
		int m = l;

		for (int i=l+1; i<=u; i++) {
			if ( the_data[ ind[i] ] [c] < the_data[t][c]) {
				m++;
				std::swap(ind[i],ind[m]); 
				}
			} // for i 
		std::swap(ind[l],ind[m]);

		if (m <= k) l = m+1;
		if (m >= k) u = m-1;
		} // while loop
	}

int kdtree2::select_on_coordinate_value(int c, MyFloat alpha, int l, int u) {
	//
	//  Move indices in ind[l..u] so that the elements in [l .. return]
	//  are <= alpha, and hence are less than the [return+1..u]
	//  elmeents, viewed across dimension 'c'.
	// 
	int lb = l, ub = u;

	while (lb < ub) {
		if (the_data[ind[lb]][c] <= alpha) {
			lb++; // good where it is.
			} else {
				std::swap(ind[lb],ind[ub]); 
				ub--;
			}
		}

	// here ub=lb
	if (the_data[ind[lb]][c] <= alpha)
		return(lb);
	else
		return(lb-1);

	}


// void kdtree2::dump_data() {
//   int upper1, upper2;

//   upper1 = N;
//   upper2 = dim;

//   printf("Rearrange=%d\n",rearrange);
//   printf("N=%d, dim=%d\n", upper1, upper2);
//   for (int i=0; i<upper1; i++) {
//     printf("the_data[%d][*]=",i); 
//     for (int j=0; j<upper2; j++)
//       printf("%f,",the_data[i][j]); 
//     printf("\n"); 
//   }
//   for (int i=0; i<upper1; i++)
//     printf("Indexes[%d]=%d\n",i,ind[i]); 
//   for (int i=0; i<upper1; i++) {
//     printf("data[%d][*]=",i); 
//     for (int j=0; j<upper2; j++)
//       printf("%f,",(*data)[i][j]); 
//     printf("\n"); 
//   }
// }



//
// search record substructure
//
// one of these is created for each search.
// this holds useful information  to be used
// during the search



static const MyFloat infinity = std::numeric_limits<float>::max();

class searchrecord {

private:
	friend class kdtree2;
	friend class kdtree2_node;

	std::vector<MyFloat>& qv; 
	int dim;
	bool rearrange;
	unsigned int nn; // , nfound;
	MyFloat ballsize;
	int centeridx, correltime;

	kdtree2_result_vector& result;  // results
	const kdtree2_array* data; 
	const std::vector<int>& ind;
	const std::vector<bool>& data_state;
	// constructor

public:
	searchrecord(std::vector<MyFloat>& qv_in, kdtree2& tree_in,
		kdtree2_result_vector& result_in) :  
	qv(qv_in),
		result(result_in),
		data(tree_in.data),
		data_state(tree_in.data_state),
		ind(tree_in.ind) 
		{
		dim = tree_in.dim;
		rearrange = tree_in.rearrange;
		ballsize = infinity; 
		nn = 0; 
		};

	};


void kdtree2::n_nearest_brute_force(std::vector<MyFloat>& qv, int nn, kdtree2_result_vector& result) {

	result.clear();

	for (int i=0; i<N; i++) {
		MyFloat dis = 0.0;
		kdtree2_result e; 
		for (int j=0; j<dim; j++) {
			dis += squared( the_data[i][j] - qv[j]);
			}
		e.dis = dis;
		e.idx = i;
		result.push_back(e);
		}
	sort(result.begin(), result.end() ); 

	}


void kdtree2::n_nearest(std::vector<MyFloat>& qv, int nn, kdtree2_result_vector& result) {
	searchrecord sr(qv,*this,result);
	//std::vector<MyFloat> vdiff(dim,0.0); 

	result.clear(); 

	sr.centeridx = -1;
	sr.correltime = 0;
	sr.nn = nn; 

	root->search(sr); 

	if (sort_results) std::sort(result.begin(), result.end());

	}
// search for n nearest to a given query vector 'qv'.
/*void kdtree2::n_density(std::vector<MyFloat>& qv, int nn, kdtree2_result_vector& result, MyFloat &den) {
	searchrecord sr(qv,*this,result);

	result.clear(); 

	sr.centeridx = -1;
	sr.correltime = 0;
	sr.nn = nn; 

	root->search(sr); 

	if (sort_results) std::sort(result.begin(), result.end());

	}
*/

void kdtree2::n_nearest_around_point(int idxin, int correltime, int nn,
									 kdtree2_result_vector& result) 
	{
	std::vector<MyFloat> qv(dim);  //  query vector

	result.clear(); 

	for (int i=0; i<dim; i++) {
		qv[i] = the_data[idxin][i]; 
		}
	// copy the query vector.

		{
		searchrecord sr(qv, *this, result);
		// construct the search record.
		sr.centeridx = idxin;
		sr.correltime = correltime;
		sr.nn = nn; 
		root->search(sr); 
		}

		if (sort_results) sort(result.begin(), result.end());

	}


void kdtree2::r_nearest(std::vector<MyFloat>& qv, MyFloat r2, kdtree2_result_vector& result) {
	// search for all within a ball of a certain radius
	searchrecord sr(qv,*this,result);
	std::vector<MyFloat> vdiff(dim,0.0); 

	result.clear(); 

	sr.centeridx = -1;
	sr.correltime = 0;
	sr.nn = 0; 
	sr.ballsize = r2; 

	root->search(sr); 

	if (sort_results) sort(result.begin(), result.end());

	}

int kdtree2::r_count(std::vector<MyFloat>& qv, MyFloat r2) {
	// search for all within a ball of a certain radius
	{
	kdtree2_result_vector result; 
	searchrecord sr(qv,*this,result);

	sr.centeridx = -1;
	sr.correltime = 0;
	sr.nn = 0; 
	sr.ballsize = r2; 

	root->search(sr); 
	return(result.size());
	}


	}

void kdtree2::r_nearest_around_point(int idxin, int correltime, MyFloat r2,
									 kdtree2_result_vector& result) 
	{
	std::vector<MyFloat> qv(dim);  //  query vector

	result.clear(); 

	for (int i=0; i<dim; i++) {
		qv[i] = the_data[idxin][i]; 
		}
	// copy the query vector.

		{
		searchrecord sr(qv, *this, result);
		// construct the search record.
		sr.centeridx = idxin;
		sr.correltime = correltime;
		sr.ballsize = r2; 
		sr.nn = 0; 
		root->search(sr); 
		}

		if (sort_results) sort(result.begin(), result.end());

	}


int kdtree2::r_count_around_point(int idxin, int correltime, MyFloat r2) 
	{
	std::vector<MyFloat> qv(dim);  //  query vector


	for (int i=0; i<dim; i++) {
		qv[i] = the_data[idxin][i]; 
		}
	// copy the query vector.

		{
		kdtree2_result_vector result; 
		searchrecord sr(qv, *this, result);
		// construct the search record.
		sr.centeridx = idxin;
		sr.correltime = correltime;
		sr.ballsize = r2; 
		sr.nn = 0; 
		root->search(sr); 
		return(result.size());
		}


	}


// 
//        KDTREE2_NODE implementation
//

// constructor
kdtree2_node::kdtree2_node(int dim) : box(dim) { 
	left = right = NULL; 
	//
	// all other construction is handled for real in the 
	// kdtree2 building operations.
	// 
	}

// destructor
kdtree2_node::~kdtree2_node() {
	if (left != NULL) delete left; 
	if (right != NULL) delete right; 
	// maxbox and minbox 
	// will be automatically deleted in their own destructors. 
	}


void kdtree2_node::search(searchrecord& sr) {
	// the core search routine.
	// This uses true distance to bounding box as the
	// criterion to search the secondary node. 
	//
	// This results in somewhat fewer searches of the secondary nodes
	// than 'search', which uses the vdiff vector,  but as this
	// takes more computational time, the overall performance may not
	// be improved in actual run time. 
	//

	if ( (left == NULL) && (right == NULL)) {
		// we are on a terminal node
		if (sr.nn == 0) {
			process_terminal_node_fixedball(sr);
			} else {
				process_terminal_node(sr);
			}
		} else {
			kdtree2_node *ncloser, *nfarther;

			MyFloat extra;
			MyFloat qval = sr.qv[cut_dim]; 
			// value of the wall boundary on the cut dimension. 
			if (qval < cut_val) {
				ncloser = left;
				nfarther = right;
				extra = cut_val_right-qval;
				} else {
					ncloser = right;
					nfarther = left;
					extra = qval-cut_val_left; 
				};

				if (ncloser != NULL) ncloser->search(sr);

				if ((nfarther != NULL) && (squared(extra) < sr.ballsize)) {
					// first cut
					if (nfarther->box_in_search_range(sr)) {
						nfarther->search(sr); 
						}      
					}
		}
	}


inline MyFloat dis_from_bnd(MyFloat x, MyFloat amin, MyFloat amax) {
	if (x > amax) {
		return(x-amax); 
		} else if (x < amin)
			return (amin-x);
		else
			return 0.0;

	}

inline bool kdtree2_node::box_in_search_range(searchrecord& sr) {
	//
	// does the bounding box, represented by minbox[*],maxbox[*]
	// have any point which is within 'sr.ballsize' to 'sr.qv'??
	//

	int dim = sr.dim;
	MyFloat dis2 =0.0; 
	MyFloat ballsize = sr.ballsize; 
	if(this->metric.set_flag){
		vector<float> vec(dim);
		for (int i=0; i<dim; i++){
			vec[i]=(float)(dis_from_bnd(sr.qv[i],box[i].lower,box[i].upper));
			}
		dis2=this->metric.doDist(vec);
		}else{
			for (int i=0; i<dim;i++) {
				dis2 += squared(dis_from_bnd(sr.qv[i],box[i].lower,box[i].upper));
				}
		if (dis2 > ballsize)
			return(false);
		}
	return(true);
	}


void kdtree2_node::process_terminal_node(searchrecord& sr) {
	int centeridx  = sr.centeridx;
	int correltime = sr.correltime;
	unsigned int nn = sr.nn; 
	int dim        = sr.dim;
	MyFloat ballsize = sr.ballsize;
	//
	bool rearrange = sr.rearrange; 
	const kdtree2_array& data = *sr.data;
	//const std::vector<bool> &data_state = *sr.data_state;

	const bool debug = false;

	if (debug) {
		printf("Processing terminal node %d, %d\n",l,u);
		cout << "Query vector = [";
		for (int i=0; i<dim; i++) cout << sr.qv[i] << ','; 
		cout << "]\n";
		cout << "nn = " << nn << '\n';
		//    check_query_in_bound(sr);
		}

	for (int i=l; i<=u;i++) {
		int indexofi;  // sr.ind[i]; 
		MyFloat dis, vol;
		bool early_exit; 
		bool is_active;
		if (rearrange) {
			early_exit = false;
			dis = 0.0;
			for (int k=0; k<dim; k++) {
				dis += squared(data[i][k] - sr.qv[k]);
				vol*=(box[k].upper-box[k].lower);
				if (dis > ballsize) {
					early_exit=true; 
					break;
					}
				}
			if(early_exit) continue; // next iteration of mainloop
			// why do we do things like this?  because if we take an early
			// exit (due to distance being too large) which is common, then
			// we need not read in the actual point index, thus saving main
			// memory bandwidth.  If the distance to point is less than the
			// ballsize, though, then we need the index.
			//
			indexofi = sr.ind[i];
			} else {
				// 
				// but if we are not using the rearranged data, then
				// we must always 
				indexofi = sr.ind[i];
				early_exit = false;
				dis = 0.0;
				vol=1.0;
				//Mahalanobis distance
				if(this->metric.set_flag){					
					vector<float> vec(dim);
					for (int k=0; k<dim; k++){
						vec[k]=(float)(data[indexofi][k]-sr.qv[k]);
						}
					dis=this->metric.doDist(vec);
					if (dis > ballsize) {
						early_exit= true; 
						break;
						}
					}else{
						for (int k=0; k<dim; k++) {
							dis += squared(data[indexofi][k] - sr.qv[k]);
							//	vol*=(box[k].upper-box[k].lower);
							if (dis > ballsize) {
								early_exit= true; 
								break;
								}
							}
					}
				if(early_exit) continue; // next iteration of mainloop
			} // end if rearrange. 

		if (centeridx > 0) {
			// we are doing decorrelation interval
			if (abs(indexofi-centeridx) < correltime) continue; // skip this point. 
			}

		// here the point must be added to the list.
		//
		// two choices for any point.  The list so far is either
		// undersized, or it is not.
		//
		if (sr.result.size() < nn) {
			kdtree2_result e;
			e.idx = indexofi;
			e.dis = dis;
			e.vol=std::abs(vol);
			e.metric=&(this->metric);
			is_active=sr.data_state[indexofi];// if point data is active in NGB searcher then put in result
			if(dis>0 && is_active)
				{
				sr.result.push_element_and_heapify(e); 
				if (debug) cout << "unilaterally pushed dis=" << dis;
				if (sr.result.size() == nn) ballsize = sr.result.max_value();
				// Set the ball radius to the largest on the list (maximum priority).
				if (debug) {
					cout << " ballsize = " << ballsize << "\n"; 
					cout << "sr.result.size() = "  << sr.result.size() << '\n';
					}
				}
			} else {
				//
				// if we get here then the current node, has a squared 
				// distance smaller
				// than the last on the list, and belongs on the list.
				// 
				kdtree2_result e;
				e.idx = indexofi;
				e.dis = dis;
				e.vol=std::abs(vol);
				e.metric=&(this->metric);
				is_active=sr.data_state[indexofi];// if point data is active in NGB searcher then put in result
				if(dis>0 && is_active)
					{
					ballsize = sr.result.replace_maxpri_elt_return_new_maxpri(e); 
					if (debug) {
						cout << "Replaced maximum dis with dis=" << dis << 
							" new ballsize =" << ballsize << '\n';
						}
					}
			}
		} // main loop
	sr.ballsize = ballsize;
	}

void kdtree2_node::process_terminal_node_fixedball(searchrecord& sr) {
	int centeridx  = sr.centeridx;
	int correltime = sr.correltime;
	int dim        = sr.dim;
	MyFloat ballsize = sr.ballsize;
	//
	bool rearrange = sr.rearrange; 
	const kdtree2_array& data = *sr.data;

	for (int i=l; i<=u;i++) {
		int indexofi = sr.ind[i];
		MyFloat dis, vol;
		bool early_exit;
		bool is_active;

		if (rearrange) {
			early_exit = false;
			dis = 0.0;
			vol=1.0;
			//Mahalanobis distance
			if(this->metric.set_flag){					
				vector<float> vec(dim);
				for (int k=0; k<dim; k++){
					vec[k]=(float)(data[indexofi][k]-sr.qv[k]);
					}
				dis=this->metric.doDist(vec);
				if (dis > ballsize) {
					early_exit= true; 
					break;
					}
				}else{
					for (int k=0; k<dim; k++) {
						dis += squared(data[i][k] - sr.qv[k]);
						//vol*=(box[k].upper-box[k].lower);
						if (dis > ballsize) {
							early_exit=true; 
							break;
							}
						}
				}
			if(early_exit) continue; // next iteration of mainloop
			// why do we do things like this?  because if we take an early
			// exit (due to distance being too large) which is common, then
			// we need not read in the actual point index, thus saving main
			// memory bandwidth.  If the distance to point is less than the
			// ballsize, though, then we need the index.
			//
			indexofi = sr.ind[i];
			} else {
				// 
				// but if we are not using the rearranged data, then
				// we must always 
				indexofi = sr.ind[i];
				early_exit = false;
				dis = 0.0;
				//Mahalanobis distance
				if(this->metric.set_flag){					
					vector<float> vec(dim);
					for (int k=0; k<dim; k++){
						vec[k]=(float)(data[indexofi][k]-sr.qv[k]);
						}
					dis=this->metric.doDist(vec);
					if (dis > ballsize) {
						early_exit= true; 
						break;
						}
					}else{
						for (int k=0; k<dim; k++) {
							dis += squared(data[indexofi][k] - sr.qv[k]);
							vol*=(box[k].upper-box[k].lower);
							if (dis > ballsize) {
								early_exit= true; 
								break;
								}
							}
					}
				if(early_exit) continue; // next iteration of mainloop
			} // end if rearrange. 

		if (centeridx > 0) {
			// we are doing decorrelation interval
			if (abs(indexofi-centeridx) < correltime) continue; // skip this point. 
			}

			{
			kdtree2_result e;
			e.idx = indexofi;
			e.dis = dis;
			e.vol=std::abs(vol);
			e.metric=&(this->metric);
			is_active=sr.data_state[indexofi];// if point data is active in NGB searcher then put in result
			if(dis>0 && is_active)
				sr.result.push_back(e);
			}

		}
	}
