#include "utils.h"
#include  <fstream>
#include <boost/lexical_cast.hpp>
using boost::lexical_cast;
using boost::bad_lexical_cast;


bool is_file_exist(const std::string &filename)
	{
	std::ifstream fin(filename.c_str(), std::ios::in);
	if(fin.is_open())
		{
		fin.close();
		return true;
		}
	fin.close();

	return false;
	}

unsigned int GetISnap(std::string str)
	{
	size_t found=str.find_last_of("_");
	unsigned int val=boost::lexical_cast<unsigned int>(str.substr(found+1));
	return val;
	}

//////////////////////////////////////
#include <cmath>
#include <vector>

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void jacobi(std::vector<std::vector<float> > &a, int n, std::vector<float> &d, std::vector<std::vector<float> > &v, int *nrot)
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c;///*,*b,*z;
	std::vector<float> b(n+1), z(n+1);
	//b=vector(1,n);
	//z=vector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0f;
		v[ip][ip]=1.0f;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0f;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0f;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0f) {
			//free_vector(z,1,n);
			//free_vector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2f*sm/(n*n);
		else
			tresh=0.0f;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0f*fabs(a[ip][iq]);
				if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
					a[ip][iq]=0.0f;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5f*h/(a[ip][iq]);
						t=1.0f/(fabs(theta)+sqrt(1.0f+theta*theta));
						if (theta < 0.0f) t = -t;
					}
					c=1.0f/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0f+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0f;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0f;
		}
	}
	std::cerr<<"Too many iterations in routine jacobi"<<std::endl;
}
#undef ROTATE

/* (C) Copr. 1986-92 Numerical Recipes Software ,2kB. */
//////////////////////////////////////