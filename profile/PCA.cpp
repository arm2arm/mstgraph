#include "PCA.h"
#include "utils.h"
CPCA::CPCA(void):verbose(false)
	{
	}

CPCA::~CPCA(void)
	{
	}

using std::vector;
template <typename T>
void print(vector<vector<T> > &Cov, int ish=1)
{
 	for(size_t i =ish;i<Cov.size()-ish;i++)
	  {
		for(size_t j =ish;j<Cov[0].size()-ish;j++)
		  {
		    cout<<Cov[i][j]<<" ";
		  }
		cout<<endl;
	  }
	cout<<endl;
}
/// Matrix tools
template <typename T> 
void Transpose(std::vector<std::vector<T> > &a, int ish=1)
	{
	for(size_t i =ish;i<a.size()-ish;i++)
	  {
		for(size_t j =ish;j<a[0].size()-ish;j++)
		  {
		  T ta=a[i][j];
		  a[i][j]=a[j][i];
		  a[j][i]=ta;
		  }
	  } 
	}

/**  Variance-covariance matrix: creation  *****************************/
void covcol(std::vector<std::vector<double> >& data,int  n,int  m, std::vector<std::vector<double> >&  symmat)
/* Create m * m covariance matrix from given n * m data matrix. */
{
int i,  j1, j2;
/* Calculate the m * m covariance matrix. */
for (j1 = 0; j1 < m; j1++)
    {
    for (j2 = j1; j2 < m; j2++)
        {
        symmat[j1][j2] = 0.0;
        for (i = 0; i < n; i++)
            {
            symmat[j1][j2] += data[i][j1] * data[i][j2];
            }
        symmat[j2][j1] = symmat[j1][j2]/(double)(n-1);
        }
    }

}
////////////////
//#include<omp.h>
void MultMat(std::vector<std::vector<double> >& m1, std::vector<std::vector<double> >& m2, std::vector<std::vector<double> >& mult)
	{
	int r1=m1.size(), c2=m2[0].size();
	int i, k, j;
	//omp_set_num_threads(2 or 4 or 6) // as depends on your loop
//#pragma omp parallel for shared(mult, r1,c2) private(j,k) 
	for(i=0;i<r1;i++)
		{
		for(j=0;j<c2;j++)
			{
			mult[i][j]=0;
			for(k=0;k<r1;k++)
				{
				mult[i][j]+=m1[i][k]*m2[k][j];
				/*mult[0][0]=m1[0][0]*m2[0][0]+m1[0][1]*m2[1][0]+m1[0][2]*m2[2][0];*/
				}
			//printf("%d\t",mult[i][j]);
			}
		//printf("\n");
		}
	}
//////////////

double CPCA::GetCovarMatrix( CMSTree &mst, int ID)
	{

	if(mst.m_MSTCatalog.size() ==0)return 0.0;	
	size_t np=mst.m_MSTCatalog[ID].size();
	// lets us get mean
	vector<double> x, y, z, eigen, eigenvalue;
	for(size_t ip=0,i;ip</**/np;ip++)
		{
			i=mst.m_MSTCatalog[ID].id[ip];
			x.push_back(mst.m_x[i]);
			y.push_back(mst.m_y[i]);
			z.push_back(mst.m_z[i]);
		}
	GetPCAXYZ(x, y, z, eigen, eigenvalue);
	//Some strange stuff is going here....
	for(size_t ii=0;ii<eigen.size();ii++)
		eigen[ii]= -eigen[ii];
	if(verbose)
		{
	cout<<"\nresult\n"<<endl;
	copy (eigen.begin(), eigen.end(),
          ostream_iterator<double>(cout," "));
		
	cout << endl;
		}
	double Phi = rad2deg<double>(atan2( (double)(eigen[3]), (double)(eigen[0])));
	cout<<"Phi = "<<Phi<<endl;
	//cout<<(double)(v[0+1][1+1])<<" "<< (double)(v[1+1][1+1])<<endl;
	//exit(0);
	return Phi;
	}

void CPCA::GetPCA(void)
	{
	size_t i;
	vector<vector<double> > Cov( 4, vector<double>(4, 0.0))/*, CovNew( 4, vector<double>(4, 0.0))*/;
	double x[] = {2.5, 0.5, 2.2, 1.9, 3.1, 2.3, 2.0, 1.0, 1.5, 1.1};
	double y[] = {2.4, 0.7, 2.9, 2.2, 3.0, 2.7, 1.6, 1.1, 1.6, 0.9};
	//double z[] = {2.4, 0.7, 2.9, 2.2, 3.0, 2.7, 1.6, 1.1, 1.6, 0.9};
	
	size_t np=sizeof(x)/sizeof(double);
	
	double m1=accumulate(x, x+np, 0.0)/(double)np;
	double m2=accumulate(y, y+np, 0.0)/(double)np;
	//double m3=accumulate(z, z+np, 0.0)/(double)np;
	std::transform( x, x+np, x,std::bind2nd( std::minus<double>(), m1) );
	std::transform( y, y+np, y,std::bind2nd( std::minus<double>(), m2) );
	//std::transform( z, z+np, y,std::bind2nd( std::minus<double>(), m3) );
		
	/*vector<vector<double> > data;
	data.resize(np);
	for(i=0;i<np;i++)data[i].resize(3,0.0);

	for(size_t i=0;i<np;i++){
		data[i][0]=x[i];
		data[i][1]=y[i];
		}*/
	//this  returns normalized to (n-1)
	//covcol( data, np, 2, CovNew);
	for(size_t ip=0;ip</**/np;ip++)
		{
			i=ip;
			Cov[0][0] += x[i]*x[i];
			Cov[0][1] += x[i]*y[i];

			Cov[1][0] += y[i]*x[i];
			Cov[1][1] += y[i]*y[i];
		}
	for(i=0;i<3;i++)
		for(size_t j=0;j<3;j++)
			Cov[i][j]/=double(np-1);

	print(Cov, 0);
	std::vector<vector<double> > a( 4, vector<double>(4, 0.0)), v( 4, vector<double>(4, 0.0));
	std::vector<double> d(4, 0.0);
	
	int n=2,nrot;
	const int ish=1;
	for(size_t i =0;i<Cov.size()-1;i++)
		for(size_t j =0;j<Cov[0].size()-1;j++)
			a[i+ish][j+ish]=Cov[j][i];//transpose here

	//print(a);
	//Transpose(a);
	print(a);
	jacobi(a,n,d,v,&nrot);
	print(a);
	print(v);
	copy (d.begin(), d.end(),
          ostream_iterator<double>(cout," "));
    cout << endl;
	double Phi = rad2deg<double>(atan2((double)(v[1+1][1+1]), (double)(v[0+1][1+1])));
	cout<<"Phi = "<<Phi<<endl;
	exit(0); 
	
	}

void CPCA::GetPCAXYZ(vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &eigenvec, vector<double> &eigenvalue)
	{
	size_t i;
	vector<vector<double> > Cov( 4, vector<double>(4, 0.0))/*, CovNew( 4, vector<double>(4, 0.0))*/;
	
	size_t np=x.size();
	
	double m1=accumulate(x.begin(), x.end(), 0.0)/(double)np;
	double m2=accumulate(y.begin(), y.end(), 0.0)/(double)np;
	double m3=accumulate(z.begin(), z.end(), 0.0)/(double)np;
	std::transform( x.begin(), x.end(), x.begin(),std::bind2nd( std::minus<double>(), m1) );
	std::transform( y.begin(), y.end(), y.begin(),std::bind2nd( std::minus<double>(), m2) );
	std::transform( z.begin(), z.end(), z.begin(),std::bind2nd( std::minus<double>(), m3) );
		
	for(size_t i=0;i<np;i++)
		{
			
			Cov[0][0] += x[i]*x[i];
			Cov[0][1] += x[i]*y[i];
			Cov[0][2] += x[i]*z[i];
	
			Cov[1][0] += y[i]*x[i];
			Cov[1][1] += y[i]*y[i];
			Cov[1][2] += y[i]*z[i];

			Cov[2][0] += z[i]*x[i];
			Cov[2][1] += z[i]*y[i];
			Cov[2][2] += z[i]*z[i];

		}
	for(i=0;i<3;i++)
		for(size_t j=0;j<3;j++)
			Cov[i][j]/=double(np-1);
if(verbose)
		
	print(Cov, 0);
	std::vector<vector<double> > a( 4, vector<double>(4, 0.0)), v( 4, vector<double>(4, 0.0));
	std::vector<double> d(4, 0.0);
	
	int n=3,nrot;
	const int ish=1;
	for(size_t i =0;i<Cov.size()-1;i++)
		for(size_t j =0;j<Cov[0].size()-1;j++)
			a[i+ish][j+ish]=Cov[i][j];//transpose here

	//print(a);
	Transpose(a);
	if(verbose)
	print(a, 0);
	jacobi(a,n,d,v,&nrot);
	if(verbose)
	print(a, 0);
	Transpose(v);
	if(verbose)
		{
	print(v, 0);
	copy (d.begin(), d.end(),
          ostream_iterator<double>(cout," "));
		}
	//cout << endl;
	//double Phi = rad2deg<double>(atan2((double)(v[1+1][1+1]), (double)(v[0+1][1+1])));
	//cout<<"Phi = "<<Phi<<endl;
	eigenvalue=d;
	for(int i=ish;i<n+ish;i++)
		for(int j=ish;j<n+ish;j++)
			eigenvec.push_back(v[i][j]);
//	cout<<"Final Result...."<<endl;
if(verbose)
	copy (eigenvec.begin(), eigenvec.end(),
          ostream_iterator<double>(cout," "));
    
//	exit(0); 
	}