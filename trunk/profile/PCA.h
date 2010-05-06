#ifndef _PCA_
#define _PCA_
#include "MSTree.h"
#include "MSTGroup.h"


class CPCA
	{
	public:
		CPCA(void);
		~CPCA(void);
		void GetCovarMatrix( CMSTree &mst, int ID=0);
    double GetPCA(void);

	private:
		//TMSTCat m_MSTCatalog;

	};



#endif
