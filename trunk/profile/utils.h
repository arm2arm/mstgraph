#ifndef _CRange_
#define _CRange_
#include <iostream>
#include <algorithm>
#include <string>
#include <cassert>
#include <vector>

using std::min;
using std::max;
using std::cout;
using std::endl;

class CRange{
public:
	float Min;
	float Max;
	float m_sum;
	float m_mean;
	unsigned int m_np;
	CRange(){Reset();};
	~CRange(){};
	void sum(float x){m_sum+=x;m_np++;};

	void getbound(float x)
	{
		Min=min(x, Min);
		Max=max(x, Max);
		sum(x);
	};
void Reset(void ){Min=1e10; Max=-1e10;m_sum=0;m_np=0;};
void print(const char *str){
	cout<<str<<" Min="<<Min<<"\tMax="<<Max<<"\t mean="<<m_sum/m_np<<endl;

};
};


std::string intToString(int i);
#define printOpenGLError() printOglError(__FILE__, __LINE__)
float getEnv(char * name, float defvalue);
int printOglError(char *file, int line);
float   stoptimer_( int *flag);
void  starttimer_( void);
float mydrand48(void);
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795f
#endif


template<typename Container>
void delete_pointers_list(Container& c) { while(!c.empty()) delete c.back(), c.pop_back(); }





template<class T>
class CKernel{
public:
        CKernel(const int dim=3):m_dim(dim){
                if((m_dim<1)||(m_dim>20))
                        {
                        std::cout<<"Specify the Normalization constant for dimensions > 20"<<endl;
                        
                        }
                AllocateTable();
                };
        ~CKernel(){};
        virtual T W(T u)=0;
protected:
        void AllocateTable()
                {
                KernelRad.resize(KERNEL_TABLE+2);//+2 allow to go in the loops  KERNEL_TABLE+1;
                Kernel.resize(KERNEL_TABLE+2);
                };
        virtual void    InitTable()=0;
protected:
        static const int KERNEL_TABLE=1000;
        std::vector<T> KernelRad,Kernel;
        int m_dim;
        std::string TypeOfKernel;
        };



template<class T>
class CEpanechikov : public CKernel<T>
        {
        public:
                CEpanechikov(int dim):CKernel<T>(dim){this->TypeOfKernel=std::string("Epanechikov kernel= Cd*(1-u)");InitTable();};
        private:
        virtual void    InitTable(){
                static double
                        fep[]={0.75000113,0.63661975,0.59683102,0.60792705,0.66492015,0.77403670,0.95242788,1.2319173,1.6674189,2.3527875,3.4499109,5.2424031,8.2360444,13.349647,22.283674,38.243824,67.384374,121.73344,225.21478,426.23651};
                T f1=(T)fep[this->m_dim-1];//Epanechikov
                cout<<"Normalization constant of Kernel type "<< this->TypeOfKernel<<": "<<f1<<endl;
                int i=0;
                for(i=0;i<=this->KERNEL_TABLE+1;i++)
                        {
                        this->KernelRad[i] = ((T)i)/(T)this->KERNEL_TABLE;
                        this->Kernel[i] = f1 *(1-this->KernelRad[i]*this->KernelRad[i]);
                        }
                this->Kernel[this->KERNEL_TABLE+1] =(T)0.0;
                };
        public:
        virtual T W(T u){
                int k;
                T wk=0.0;
                if(u<(T)1.0)
                        {
                        k = (int)(u*this->KERNEL_TABLE);
                        wk =( this->Kernel[k]  + (this->Kernel[k+1]-this->Kernel[k])*(u-this->KernelRad[k])*this->KERNEL_TABLE);
                        //rhoxyz+=wk*Part[pqStartA[i].p].Mass;
                        }
                return wk;
                };

        };



bool is_file_exist(const std::string &filename);
unsigned int GetISnap(std::string str);
/////////////////////
template <class T>
class ProxyLess
	{
	T& that;
	public:
		ProxyLess(T &f) : that(f){}
		bool operator()(unsigned int leftID, unsigned int rightID) const 
			{
			return that[leftID]<that[rightID];
			}
	};
///////////////////////
#include <boost/timer.hpp>
#include <boost/config.hpp>
#include <ctime>
#include <cstring>

#include <boost/limits.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"

class scoped_timer {
	boost::posix_time::ptime start_;
	std::string m_text;
public:    
	inline  void   SetText(std::string text){m_text=text;};
	scoped_timer(std::string text) 
		: start_(boost::posix_time::microsec_clock::universal_time()),m_text(text)
		{
		//m_start= getticks();
		}
	~scoped_timer() {

		boost::posix_time::ptime stop( boost::posix_time::microsec_clock::universal_time() );
//		m_end=getticks();
//		double clock_cycles=elapsed(m_end, m_start);
		std::cout<<std::setprecision(2)<<std::fixed;
		std::cout <<"# "<<m_text<< " done in " << ( stop - start_).total_milliseconds() << " milli seconds "<<std::endl;
		std::cout<<std::setprecision(16)<<std::fixed;
		}
protected:
	
	};
///////////////////////////////////////////
size_t count_nonblanks(std::string &str);
///////////////// USED by REGEX ///////////
template<class In,class Out,class Pred, class Op>
Out transform_if(In first,In last,Out res,Pred p,Op op)
	{
	while(first!=last){
		if (p(*first))
			*res = op(*first);
		++first;++res;
		}
	return res;
	}

////////////////////////
template <class T >
T rad2deg(T v){return v*(T)(180/M_PI);};
template <class T >
T deg2rad(T v){return v*(T)(M_PI/180.0);};

///////////////////////
#define SafeFree(a) if(a!=NULL)delete a;


//void jacobi(std::vector<std::vector<float> > &a, int n, std::vector<float>  &d, std::vector<std::vector<float> > &v, int *nrot);
//////////////////////////////////////
#include <cmath>
#include <vector>

void jacobi(std::vector<std::vector<double> > &a, int n, std::vector<double> &d, std::vector<std::vector<double> > &v, int *nrot, int ish=1);
void eigsrt(std::vector<double>  &d, std::vector<std::vector<double> > &v, int n);

//////////////////////////////////////
using std::vector;
using std::cout;
using std::endl;
template <typename T>
void dump_xyz(std::vector<T>&x,std::vector<T> &y,vector<T> &z, std::string filename="dum_points.log")
	{
	cout<<"dumping";
	std::string fname=filename;
	std::ofstream of(fname.c_str());
	if(of.is_open())
		for(size_t ip=0;ip<x.size();ip++)
			{
			if((ip%10000)==0)
				cout<<".";
			cout.flush();
			of<<std::setw(12)<<std::setprecision(7)<<x[ip]<<" "<<y[ip]<<" "<<z[ip]<<endl;

			}
		of.close();
		cout<<"done!"<<endl;
	}


template<class T> 
std::vector<T> smooth(std::vector<T> &v)
	{

	if(v.size()<5)return v;
	int np=v.size();
	std::vector<T> vsm=v;

	for(int i=2;i<np-2;i++)
		vsm[i]=(v[i-2]+2*v[i-1]+3*v[i]+2*v[i+1]+v[i+2])/9.0f;
	vsm[0]=(vsm[0]+vsm[1]+vsm[2])/3.0f;
	vsm[1]=(vsm[1]+vsm[2]+vsm[3])/3.0f;
	vsm[np-2]=(vsm[np-2]+vsm[np-3]+vsm[np-4])/3.0f;
	vsm[np-1]=(vsm[np-1]+vsm[np-2]+vsm[np-3])/3.0f;
	return vsm;
	}
template<class T>
std::vector< std::vector<T> > GetAB(vector<T> &R,vector<T> &x,vector<T> &y, T profRMAX=10.0, T xStep=0.1)
	{
	std::vector<unsigned int> idxR(R.size(), 0);
	for(size_t ii=0;ii<R.size(); ii++)
		{
		idxR[ii]=ii;
		}

	std::sort(idxR.begin(), idxR.end(), ProxyLess< vector<float> >(R)); //sorting by Radius;

	/////////////////////////////////////////////
	//GetAB
	std::vector< vector<T> > result;
	int fm, ii;
	int iR20 = std::count_if(R.begin(), R.end(),
		std::bind2nd(std::less_equal<T>(), profRMAX));

	int nshell = static_cast<int> (profRMAX / xStep) + 1;
	for (fm = 2; fm < 9; fm += 2) {
		vector<T> al, al0, al2, bl2, fR;
		al.resize(nshell);
		al2.resize(nshell);
		bl2.resize(nshell);
		al0.resize(nshell);
		fR.resize(nshell);
		std::fill(al.begin(), al.end(), (T)0.0);
		std::fill(al0.begin(), al0.end(),(T) 0.0);
		std::fill(al2.begin(), al2.end(), (T)0.0);
		std::fill(bl2.begin(), bl2.end(), (T)0.0);
		std::fill(fR.begin(), fR.end(), (T)0.0);

		int ishell = 0;
		T theta = 0.0;
		T massp = 1.0;


		for (int i = 0; i < iR20; i++) {
			ii=idxR[i];
			ishell = static_cast<int> (R[ii] / xStep);
			theta = atan2(y[ii], x[ii]);
			al2[ishell] += massp * cos(fm * theta);
			bl2[ishell] += massp * sin(fm * theta);
			al0[ishell] += massp;
			}

		for (int i = 0; i < nshell; i++) {
			al[i] = sqrt(al2[i] * al2[i] +
				bl2[i] * bl2[i]);
			if (al0[i] > 0.0f)
				al[i] /= al0[i];
			fR[i] = xStep*i;
			}
		double total = std::accumulate(al.begin(), al.end(), 0.0);
		//string fmFile = fname + string("_Afm") + ToString(fm) + '_' + ToString(TYPE);
		if(fm == 2)	result.push_back(fR);
		result.push_back(smooth<float>(al));
		//result.push_back(al);
		//result.push_back(al2);
		//result.push_back(bl2);
		}
	return result;
	}



#endif

