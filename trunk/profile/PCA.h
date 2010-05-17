#ifndef _PCA_
#define _PCA_
#include "MSTree.h"
#include "MSTGroup.h"


class CPCA
	{
	public:
		CPCA(void);
		~CPCA(void);
		double  GetCovarMatrix( CMSTree &mst, int ID=0);
		void GetPCA(void){};
		void eigen (vector<vector<double> > &a, vector<double> &eigen, vector<vector<double> > &eigenvec);
void GetPCAXYZ(vector<double> &x, vector<double> &y, vector<double> &z);
inline double Phi(){return m_Phi;};
	private:
		vector<double>  eigenvals;
		std::vector<std::vector<double> > eigenvec;
		double m_Phi;
bool verbose;
	};



#endif
