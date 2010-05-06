#include "PCA.h"
#include "utils.h"
CPCA::CPCA(void)
	{
	}

CPCA::~CPCA(void)
	{
	}

using std::vector;
void print(vector<vector<float> > &Cov)
{
 	for(size_t i =1;i<Cov.size()-1;i++)
	  {
		for(size_t j =1;j<Cov[0].size()-1;j++)
		  {
		    cout<<Cov[i][j]<<" ";
		  }
		cout<<endl;
	  }
	cout<<endl;
}

double CPCA::GetCovarMatrix( CMSTree &mst, int ID)
	{
	size_t i;
	vector<vector<float> > Cov( 4, vector<float>(4, 0.0f));
	for(size_t ip=0;i<mst.m_MSTCatalog[ID].size();ip++)
		{
		i=mst.m_MSTCatalog[ID].id[ip];
		//grp.
		Cov[0][0] += mst.m_x[i]*mst.m_x[i];
		Cov[0][1] += mst.m_x[i]*mst.m_y[i];
		//Cov[0][1] += mst.m_x[i]*mst.m_y[i];
	
		Cov[1][0] += mst.m_y[i]*mst.m_x[i];
		Cov[1][1] += mst.m_y[i]*mst.m_y[i];
		
		//CovVec[0].x+=G1[i].x*G1[i].x;
		//CovVec[0].y+=G1[i].x*G1[i].y;

		}
	for(size_t i =0;i<Cov.size();i++)
		for(size_t j =0;j<Cov[0].size();j++)
			Cov[i][j]=Cov[i][j]/(float)mst.m_x.size();

	std::vector<vector<float> > a( 4, vector<float>(4, 0.0f)), v( 4, vector<float>(4, 0.0f));
	std::vector<float> d(4, 0.0f);

	int n=2,nrot;
	a[1][1] = Cov[0][0];
	a[2][2] = Cov[1][1];
	a[1][2] = Cov[0][1];
	a[2][1] = Cov[1][0];

	//	a[1][1]=0.037727899;   a[1][2]= 0.021687787;
	//	a[2][1]= 0.021687787;  a[2][2]=  0.24034810;
	/*a(1,3) = a13
	a(3,1) = a13
	a(2,3) = a23
	a(3,2) = a23
	*/
	//	print(a);
	jacobi(a,n,d,v,&nrot);
	//	print(a);
	//	print(v);
	double Phi = rad2deg<double>(atan2((double)(v[1+1][1+1]), (double)(v[0+1][1+1])));
	//cout<<"Phi = "<<Phi<<endl;
	return Phi;
	}

void CPCA::GetPCA(void)
	{
/*     SUBROUTINE get_axes(a11,a22,a33,a12,a13,a23,root1,root2,root3)
    real, INTENT(INOUT)::a11,a22,a33,a12,a13,a23
    real, INTENT(OUT)::root1,root2,root3
    integer n,idx(3)
    real a(3,3),d(3),v(3,3),tmp(3)

    n = 3
*/

/*	int n=2,nrot;
    a(1,1) = a11
    a(2,2) = a22
    a(3,3) = a33
    a(1,2) = a12
    a(2,1) = a12
    a(1,3) = a13
    a(3,1) = a13
    a(2,3) = a23
    a(3,2) = a23



	std::vector<vector<float> > a( 4, vector<float>(4, 0.0f)), v( 4, vector<float>(4, 0.0f));
	std::vector<float> d(4, 0.0f);
	
	jacobi(a,n,d,v,&nrot);*/
//"Before eigen vectors"
/*    call jacobi(a,n,n,d,v,nrot)
    root1=d(1);root2=d(2);root3=d(3)

!!!!!!!! Sort them !!!!!!!!!!!!!
    tmp(1)=max(root1, max(root2, root3))
    tmp(3)=min(root1, min(root2, root3))
    if( (root1.ne.tmp(1)) .and.(root1.ne.tmp(3))) tmp(2)=root1
    if( (root2.ne.tmp(1)) .and.(root2.ne.tmp(3))) tmp(2)=root2
    if( (root3.ne.tmp(1)) .and.(root3.ne.tmp(3))) tmp(2)=root3

    root1=tmp(1)
    root2=tmp(2)
    root3=tmp(3)

!!!!!!!!!!!END SORTING !!!!!!!!!
    a11 = v(1,1)
    a22 = v(2,2)
    a33 = v(3,3)
    a12 = v(1,2)
    a13 = v(1,3)
    a23 = v(2,3)

*/

	
	}
