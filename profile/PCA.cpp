#include "PCA.h"
#include "utils.h"
#include  "Render.h"
#include <vector>
using std::vector;
//////////////////////////
#define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) )
#define SQR(a)  ( (a)*(a) )
double pythag(double a, double b)
	{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
	}
//////////////////////////
void tqli(vector<double> &d, vector<double> &e, int n, vector<vector<double> > &z)
	{
	//float pythag(float a, float b);
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i];
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
				}
			if (m != l) {
				if (iter++ == 30) cerr<<"Too many iterations in tqli"<<endl;
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
						}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=1;k<=n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
						}
					}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
				}
			} while (m != l);
		}
	}

//////////////////////////



void tred2(vector<vector<double> >&a, int n,vector<double> &d, vector<double>&e)
/* Householder reduction of matrix a to tridiagonal form.
Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
Springer-Verlag, 1976, pp. 489-494.
W H Press et al., Numerical Recipes in C, Cambridge U P,
1988, pp. 373-374.  */
	{
	int l, k, j, i;
	double  scale, hh, h, g, f;

	for (i = n; i >= 2; i--)
		{
		l = i - 1;
		h = scale = 0.0;
		if (l > 1)
			{
			for (k = 1; k <= l; k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i] = a[i][l];
			else
				{
				for (k = 1; k <= l; k++)
					{
					a[i][k] /= scale;
					h += a[i][k] * a[i][k];
					}
				f = a[i][l];
				g = f>0 ? -sqrt(h) : sqrt(h);
				e[i] = scale * g;
				h -= f * g;
				a[i][l] = f - g;
				f = 0.0;
				for (j = 1; j <= l; j++)
					{
					a[j][i] = a[i][j]/h;
					g = 0.0;
					for (k = 1; k <= j; k++)
						g += a[j][k] * a[i][k];
					for (k = j+1; k <= l; k++)
						g += a[k][j] * a[i][k];
					e[j] = g / h;
					f += e[j] * a[i][j];
					}
				hh = f / (h + h);
				for (j = 1; j <= l; j++)
					{
					f = a[i][j];
					e[j] = g = e[j] - hh * f;
					for (k = 1; k <= j; k++)
						a[j][k] -= (f * e[k] + g * a[i][k]);
					}
				}
			}
		else
			e[i] = a[i][l];
		d[i] = h;
		}
	d[1] = 0.0;
	e[1] = 0.0;
	for (i = 1; i <= n; i++)
		{
		l = i - 1;
		if (d[i])
			{
			for (j = 1; j <= l; j++)
				{
				g = 0.0;
				for (k = 1; k <= l; k++)
					g += a[i][k] * a[k][j];
				for (k = 1; k <= l; k++)
					a[k][j] -= g * a[k][i];
				}
			}
		d[i] = a[i][i];
		a[i][i] = 1.0;
		for (j = 1; j <= l; j++)
			a[j][i] = a[i][j] = 0.0;
		}
	}

/**  Tridiagonal QL algorithm -- Implicit  **********************/

void tqli1(vector<double> &d, vector<double> &e,int n, vector<vector<double> > &z)
	{
	int m, l, iter, i, k;
	double s, r, p, g, f, dd, c, b;

	for (i = 2; i <= n; i++)
		e[i-1] = e[i];
	e[n] = 0.0;
	for (l = 1; l <= n; l++)
		{
		iter = 0;
		do
			{
			for (m = l; m <= n-1; m++)
				{
				dd = fabs(d[m]) + fabs(d[m+1]);
				if (fabs(e[m]) + dd == dd) break;
				}
			if (m != l)
				{
				if (iter++ == 30) cerr<<"No convergence in TLQI."<<endl;
				g = (d[l+1] - d[l]) / (2.0 * e[l]);
				r = sqrt((g * g) + 1.0);
				g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i = m-1; i >= l; i--)
					{
					f = s * e[i];
					b = c * e[i];
					if (fabs(f) >= fabs(g))
						{
						c = g / f;
						r = sqrt((c * c) + 1.0);
						e[i+1] = f * r;
						c *= (s = 1.0/r);
						}
					else
						{
						s = f / g;
						r = sqrt((s * s) + 1.0);
						e[i+1] = g * r;
						s *= (c = 1.0/r);
						}
					g = d[i+1] - p;
					r = (d[i] - g) * s + 2.0 * c * b;
					p = s * r;
					d[i+1] = g + p;
					g = c * r - b;
					for (k = 1; k <= n; k++)
						{
						f = z[k][i+1];
						z[k][i+1] = s * z[k][i] + c * f;
						z[k][i] = c * z[k][i] - s * f;
						}
					}
				d[l] = d[l] - p;
				e[l] = g;
				e[m] = 0.0;
				}
			}  while (m != l);
		}
	}
///////////////////////////////
#undef SIGN

CPCA::CPCA(void):verbose(true)
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
		//cout<<" eigenvector["<<i<<"]=";
		for(size_t j =ish;j<Cov[0].size()-ish;j++)
			{
			cout<<std::setw(10)<<std::setprecision(6)<<Cov[i][j]<<"  ";
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
	GetPCAXYZ(x, y, z);
	return Phi();
	}
void eigen (vector<double> &data);

void CPCA::GetPCAXYZ(vector<double> &x, vector<double> &y, vector<double> &z)
	{
	size_t i;
	vector<vector<double> > Cov( 4, vector<double>(4, 0.0))/*, CovNew( 4, vector<double>(4, 0.0))*/;

	size_t np=x.size();

	double m1=accumulate(x.begin(), x.end(), 0.0)/(double)np;
	double m2=accumulate(y.begin(), y.end(), 0.0)/(double)np;
	double m3=accumulate(z.begin(), z.end(), 0.0)/(double)np;
	cout<<m1<<" "<<m2<<" "<<m3<<" "<<endl;

	std::transform( x.begin(), x.end(), x.begin(),std::bind2nd( std::minus<double>(), m1) );
	std::transform( y.begin(), y.end(), y.begin(),std::bind2nd( std::minus<double>(), m2) );
	std::transform( z.begin(), z.end(), z.begin(),std::bind2nd( std::minus<double>(), m3) );
	vector<double> dtest;
	for(size_t i=0;i<np;i++)
		{
					
		Cov[0][0] += x[i]*x[i];
		Cov[1][0] += x[i]*y[i];
		Cov[2][0] += x[i]*z[i];

		Cov[0][1] += y[i]*x[i];
		Cov[1][1] += y[i]*y[i];
		Cov[2][1] += y[i]*z[i];

		Cov[0][2] += z[i]*x[i];
		Cov[1][2] += z[i]*y[i];
		Cov[2][2] += z[i]*z[i];

		}
	for(i=0;i<3;i++)
		for(size_t j=0;j<3;j++)
			{
			Cov[i][j]/=double(np-1);
			}
		if(verbose)		
			print(Cov, 0);

		std::vector<vector<double> > a( 4, vector<double>(4, 0.0)), v( 4, vector<double>(4, 0.0));
		std::vector<double> d(4, 0.0);

		size_t n=3;
		const int ish=1;
		for(size_t i =0;i<Cov.size()-1;i++)
			for(size_t j =0;j<Cov[0].size()-1;j++)
				{
				a[i+ish][j+ish]=Cov[j][i];// Transpose
				}

			if(verbose)
				print(a,0);
			////
			eigenvals.resize(n+1);
			/*eigenvec.resize(n+1, 0.0);
			for(size_t i=0;i<n+1;i++)
				eigenvec[i].resize(0.0);
*/
			eigenvec=a;
			eigen (a, eigenvals, eigenvec);

			cout<<"Vectors"<<endl;
			print(eigenvec, 0);
			cout<<"Values"<<endl;
			if(verbose)
				copy (eigenvals.begin(), eigenvals.end(),
				ostream_iterator<double>(cout," "));
			m_Phi=rad2deg<double>( atan2(eigenvec[1][2], eigenvec[1][1]) );
			cout<<"\nPhi="<<m_Phi<<"\n"<<endl;
			//exit(0); 
	}


void CPCA::eigen (vector<vector<double> > &a, vector<double> &eigenvalue, vector<vector<double> > &eigenvec)
	{
	size_t n=a.size()-1;
	vector<double> interm(n+1);
	vector<vector<double> > symmat=a;
	tred2(symmat, n, eigenvalue, interm);  /* Triangular decomposition */
	tqli(eigenvalue, interm, n, symmat);   /* Reduction of sym. trid. matrix */
	eigsrt(eigenvalue,symmat,n);// Let us sort them
	
	for(size_t i=1;i<=3;i++)
		{
		for(size_t j=1;j<=3;j++)
			{
			eigenvec[j][i]=symmat[i][j];//trtanspose to get right order
			}
		if(verbose)cout<<endl;
		}
	//////////////////////////

	}    



