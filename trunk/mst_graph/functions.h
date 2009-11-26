#ifndef _MYFUNCTIONS_
#define _MYFUNCTIONS_
// SPH KERNEL///////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>
#include <boost/timer.hpp>
#include <boost/config.hpp>
#include <ctime>
#include "cycle.h"
#include <boost/limits.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"
//using namespace boost;
template<class T>
class CKernel{
public:
	CKernel(const int dim=3):m_dim(dim){
		if((m_dim<1)||(m_dim>20))
			{
			cout<<"Specify the Normalization constant for dimensions > 20"<<endl;
			assert(0);
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
	virtual void	InitTable()=0;
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
		CEpanechikov(int dim):CKernel<T>(dim){TypeOfKernel="Epanechikov kernel= Cd*(1-u)";InitTable();};
	private:
	virtual void	InitTable(){
		static double
			fep[]={0.75000113,0.63661975,0.59683102,0.60792705,0.66492015,0.77403670,0.95242788,1.2319173,1.6674189,2.3527875,3.4499109,5.2424031,8.2360444,13.349647,22.283674,38.243824,67.384374,121.73344,225.21478,426.23651};
		T f1=(T)fep[m_dim-1];//Epanechikov
		cout<<"Normalization constant of Kernel type "<< TypeOfKernel<<": "<<f1<<endl;
		int i=0;
		for(i=0;i<=KERNEL_TABLE+1;i++)
			{
			KernelRad[i] = ((T)i)/(T)KERNEL_TABLE;
			Kernel[i] = f1 *(1-KernelRad[i]*KernelRad[i]);
			}
		Kernel[KERNEL_TABLE+1] =(T)0.0;
		};
	virtual T W(T u){
		int k;
		T wk=0.0;
		if(u<(T)1.0)
			{
			k = (int)(u*KERNEL_TABLE);
			wk =( Kernel[k]  + (Kernel[k+1]-Kernel[k])*(u-KernelRad[k])*KERNEL_TABLE);      
			//rhoxyz+=wk*Part[pqStartA[i].p].Mass;
			}
		return wk;
		};

	};

template<class TF>
TF Wsph(TF rr, TF h)
	{
	TF retval=0;
	TF u=rr/h;
	if(0<=u&&u<=1 )
		retval=1.0f-3.0f/2.0f*u*u+3.0f/4.0f*u*u*u;
	else
		if(1<u&&u<=2)
			retval=1.0f/4.0f*(2.0f-u)*(2.0f-u)*(2.0f-u);
		else
			return 0;
	return retval/(3.1456f*h*h*h);
	}
///////////////////////////////

template <typename T>
std::string toString(const T &thing) {
	std::ostringstream os;
	os << thing;
	return os.str();
	};

class scoped_timer {
	boost::posix_time::ptime start_;
	std::string m_text;
public:    
	inline  void   SetText(std::string text){m_text=text;};
	scoped_timer(std::string text) 
		: start_(boost::posix_time::microsec_clock::universal_time()),m_text(text)
		{
		m_start= getticks();
		}
	~scoped_timer() {

		boost::posix_time::ptime stop( boost::posix_time::microsec_clock::universal_time() );
		m_end=getticks();
		double clock_cycles=elapsed(m_end, m_start);
		std::cout<<std::setprecision(2)<<std::fixed;
		std::cout <<" "<<m_text<< " done in " << ( stop - start_).total_milliseconds() << " milli seconds or "<<clock_cycles<<" CPU cycles "<<std::endl;
		std::cout<<std::setprecision(16)<<std::fixed;
		}
protected:
	ticks m_start;
	ticks m_end;
	};

#endif