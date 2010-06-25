#ifndef _PCA_
#define _PCA_
#include "MSTree.h"
#include "MSTGroup.h"

using std::copy;
using std::ofstream;
using std::ostream_iterator;
class CPCA
	{
	public:
		CPCA(void);
		~CPCA(void);
		double  GetCovarMatrix( CMSTree &mst, int ID=0);
		void GetPCA(void){};
		void eigen (vector<vector<double> > &a, vector<double> &eigen, vector<vector<double> > &eigenvec);
void GetPCAXYZ(vector<double> &x, vector<double> &y, vector<double> &z);
void GetPCAXYZVEL(vector<double> &x, vector<double> &y, vector<double> &z,vector<double> &vx, vector<double> &vy, vector<double> &vz);
inline double Phi(){return m_Phi;};
 void dump(std::string file="eigendata.log")
   {
       std::string fname=file;
       ofstream of(fname.c_str());
       if(of.is_open())
	 {
	   of<<std::setw(10)<<std::setprecision(7);
	   copy (eigenvals.begin()+1, eigenvals.end(),
		   ostream_iterator<double>(of,"\t"));
	   copy (eigenvec[1].begin()+1, eigenvec[1].end(),
		   ostream_iterator<double>(of,"\t"));
	   copy (eigenvec[2].begin()+1, eigenvec[2].end(),
		 ostream_iterator<double>(of,"\t"));
	   copy (eigenvec[3].begin()+1, eigenvec[3].end(),
		 ostream_iterator<double>(of,"\t"));
	   
	   of<<endl;
	 }
       of.close();
   } ;
	public:
		std::vector<double>  eigenvals;
		std::vector<std::vector<double> > eigenvec;
		double m_Phi;
bool verbose;
	};



#endif
