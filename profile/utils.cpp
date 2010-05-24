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

size_t count_nonblanks(std::string &line)
	{
	size_t numSpaces = static_cast<std::size_t>(std::count(line.begin(), line.end(), ' '));
	numSpaces += static_cast<std::size_t>(std::count(line.begin(), line.end(), '\t'));
	return numSpaces;
	}


///////////////////////////////////////////

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void jacobi(std::vector<std::vector<double> > &a, int n, std::vector<double> &d, std::vector<std::vector<double> > &v, int *nrot, int ish)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c;///*,*b,*z;
	std::vector<double> b(n+ish), z(n+ish);
	for (ip=ish;ip<=(n-1+ish);ip++) {
		for (iq=ish;iq<=(n-1+ish);iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=ish;ip<=(n-1+ish);ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=ish;ip<=(n-2+ish);ip++) {
			for (iq=ip+ish;iq<=(n-1+ish);iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			//free_vector(z,1,n);
			//free_vector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=ish;ip<=(n-2+ish);ip++) {
			for (iq=ip+ish;iq<=(n-1+ish);iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
					&& (fabs(d[iq])+g) == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=ish;j<=(ip-2+ish);j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+ish;j<=(iq-2+ish);j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+ish;j<=(n-1+ish);j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=ish;j<=(n-1+ish);j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=ish;ip<=(n-1+ish);ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	std::cerr<<"Too many iterations in routine jacobi"<<std::endl;
}
#undef ROTATE
void eigsrt(std::vector<double>  &d, std::vector<std::vector<double> > &v, int n)
{
        int k,j,i;
        double p;

        for (i=1;i<n;i++) {
                p=d[k=i];
                for (j=i+1;j<=n;j++)
                        if (d[j] >= p) p=d[k=j];
                if (k != i) {
                        d[k]=d[i];
                        d[i]=p;
                        for (j=1;j<=n;j++) {
                                p=v[j][i];
                                v[j][i]=v[j][k];
                                v[j][k]=p;
                        }   
                }   
        }   
}

/* (C) Copr. 1986-92 Numerical Recipes Software ,2kB. */
//////////////////////////////////////
