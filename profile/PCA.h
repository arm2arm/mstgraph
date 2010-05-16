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
void GetPCAXYZ(vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &eigen, vector<double> &eigenvalue);
	private:
		//TMSTCat m_MSTCatalog;
bool verbose;
	};



#endif
