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
		void GetPCA(void);

	private:
		//TMSTCat m_MSTCatalog;

	};



#endif
