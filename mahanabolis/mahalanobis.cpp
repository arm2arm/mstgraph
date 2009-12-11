// Mahalanobis.cpp : Defines the entry point for the console application.
//
//http://people.revoledu.com/kardi/tutorial/Similarity/MahalanobisDistance.html
#include "stdafx.h"
  
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric::ublas;
#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include <cassert>

#include "matrix_tools.h"
#include "metric_factory.hpp"
/////////////////
/*
   Transpose of a square matrix, do it in place
*/
void Transpose(double **a,int n)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}
///////////////
class TData;
typedef std::vector<TData> TVecData;
///////////////
class TData{
public:
	TData(float x_=0.0f, float y_=0.0f):x(x_),y(y_){};	
	float x,y;
	TData operator+(TData op2){
		TData temp; 
		temp.x = x + op2.x;  
		temp.y = y + op2.y;    
		return temp; 
		};
	TData operator-(TData op2){
		TData temp; 
		temp.x = x - op2.x;  
		temp.y = y - op2.y;    
		return temp; 
		};
	TData operator+(float op2){
		TData temp; 
		temp.x = x + op2;  
		temp.y = y + op2;    
		return temp; 
		};
	TData operator-(float op2){
		TData temp; 
		temp.x = x - op2;  
		temp.y = y - op2;    
		return temp; 
		};
	TData operator/(float op2){
		TData temp; 
		temp.x = x/op2;  
		temp.y = y/op2;    
		return temp; 
		};
	TData operator*(float op2){
		TData temp; 
		temp.x = x*op2;  
		temp.y = y*op2;    
		return temp; 
		};

	};

template<typename MatrixT,typename value_type>
value_type Mean(const MatrixT& matrix)
	{
	value_type sum=std::accumulate(matrix.begin(), matrix.end(), value_type());		
	value_type mean=sum/(float)matrix.size();		
	return mean;
	};
void SwapData(float *a, float *b)
	{
		float c=(*a);
		(*a)=(*b);
		(*b)= c;
	}
template<typename MatrixT>
MatrixT Transpose(const MatrixT& matrix)
	{
	MatrixT matrixT=matrix;
	size_t nvec=matrix.size();
	for(size_t i=0;i<nvec;i++)
		{
		std::swap(matrixT[i].x,matrixT[nvec-1-i].y);
		}
	return matrixT;
	};

struct adder{
	TData m_toad;
	adder(TData toad):m_toad(toad){};
  TData operator() (TData indat) {return indat=indat+m_toad;}
};
//////////////////////////////////////////////////
TVecData GetCovarMatrix(TVecData G1)
	{
	TVecData CovVec(2);
	TVecData G1T=Transpose<TVecData>(G1);

	for(size_t i=0;i<G1.size();i++)
		{
		CovVec[0].x+=G1[i].x*G1[i].x;
		CovVec[0].y+=G1[i].x*G1[i].y;

		CovVec[1].x+=G1[i].y*G1[i].x;
		CovVec[1].y+=G1[i].y*G1[i].y;

		}
	CovVec[0]=CovVec[0]/(float)G1.size();
	CovVec[1]=CovVec[1]/(float)G1.size();
	return CovVec;
	}
TVecData operator+( TVecData& a,  TVecData& b)
	{
	if(a.size() != b.size())assert(0);
	TVecData result(a);
	for(size_t i=0;i<a.size();i++)
		result[i]=a[i]+b[i];
	return result;
	}
TVecData operator*( TVecData& a,  float b)
	{
	TVecData result(a);
	for(size_t i=0;i<a.size();i++)
		result[i]=a[i]*b;
	return result;
	}
TVecData sqrt( TVecData& a)
	{
		TVecData result(a);
		for(size_t i=0;i<a.size();i++)
			result[i]=TData(sqrt(a[i].x),sqrt(a[i].y)) ;
	return result;
	}

TVecData Inverse2x2( TVecData& a)
	{
	TVecData result(a);
	result[0].x=  a[1].y; result[0].y= -a[0].y;
	result[1].x= -a[1].x; result[1].y= -a[0].x;
	float det=a[0].x*a[1].y-a[0].y*a[1].x;
	return result*(1.0f/det);
	};

template <class T>
void GetMeanByColumns(const boost::numeric::ublas::matrix<T> &m, boost::numeric::ublas::vector<T>& vout)
	{
	// loop over columns getting each row 
	vout.resize(m.size1());
	for (size_t i = 0; i < m.size1(); ++i) 
		{ 
		vout(i)=0.0;		
		for (size_t j = 0; j < m.size2(); ++ j)
			vout(i) += m(i,j);		
		vout(i)/=(double)m.size2();
		}
	
	}
	//////////////////////////////////////////////////
template <class T>
void GetMeanByColumns_bad(const boost::numeric::ublas::matrix<T> &m, boost::numeric::ublas::vector<T>& vout)
	{
	typedef boost::numeric::ublas::matrix<T>::const_iterator1 it1;
	typedef boost::numeric::ublas::matrix<T>::const_iterator2 it2;
	size_t i=0;	
	vout.resize(m.size1());
	for (it1 i1 = m.begin1(); i1 != m.end1(); i1++)
		{	
		vout(i)=0.0;
		for (it2 i2 = m.begin2(); i2 != m.end2(); ++i2)
			{
			vout(i) += (*i2);
			}
		vout(i) /=(double)m.size2();
		i++;
		}
	}
template <class T>
boost::numeric::ublas::matrix<T> SubstractMeans( boost::numeric::ublas::matrix<T> &m, boost::numeric::ublas::vector<T> &vMeans)
	{
		ublas::matrix<T> mout(m.size1(), m.size2());
		for(size_t i=0;i<m.size1();i++)
			for(size_t j=0;j<m.size2();j++)
			{
				mout(i,j)=m(i, j)-vMeans(i);				
			}
		//cout<<mout<<endl;
		return mout;
	}
template <class T>
T pow2(T indat) {return indat*indat;};

// Mahalanobis distance!!!
int _tmain(int argc, _TCHAR* argv[])
	{

	TVecData G1,G2, Go1, Go2;

	G1.push_back(TData(2,2));
	G1.push_back(TData(2,5));
	G1.push_back(TData(6,5));
	G1.push_back(TData(7,3));
	G1.push_back(TData(4,7));
	G1.push_back(TData(6,4));
	G1.push_back(TData(5,3));
	G1.push_back(TData(4,6));
	G1.push_back(TData(2,5));
	G1.push_back(TData(1,3));
	Go1=G1;
	G2.push_back(TData(6,5));
	G2.push_back(TData(7,4));
	G2.push_back(TData(8,7));
	G2.push_back(TData(5,6));
	G2.push_back(TData(5,4));
	Go2=G2;
/////////////////////////
	TData meanG1=Mean<TVecData,TData>(G1);
	TData meanG2=Mean<TVecData,TData>(G2);
	
/////////////////////////
	std::transform( G1.begin(), G1.end(), G1.begin(), adder(TData()-meanG1) );// just pass the negated Mean
	std::transform( G2.begin(), G2.end(), G2.begin(), adder(TData()-meanG2) );// just pass the negated Mean
/////////////////////////	
	
	TVecData G1Cov=GetCovarMatrix(G1);
	TVecData G2Cov=GetCovarMatrix(G2);
	TVecData pooledCovar(G1Cov);
	float n1=(float)G1.size(), n2=(float)G2.size();
	pooledCovar=(G1Cov*(n1)+G2Cov*(n2))*(1.0f/(n1+n2));
	TVecData InvCovar=Inverse2x2(pooledCovar);
	TData mean_diff=meanG1-meanG2;

	using namespace boost::numeric::ublas;
    matrix<double> m (2, 2);
	m(0,0)=InvCovar[0].x;
	m(0,1)=InvCovar[0].y;
	m(1,0)=InvCovar[1].x;
	m(1,1)=InvCovar[1].y;
	matrix<double> meand(1, 2);
	meand(0,0)=mean_diff.x;
	meand(0,1)=mean_diff.y;
	
	
	///get covariance matrix using UBLAS
	matrix<double> uG1(2,Go1.size()),uG2(2,Go2.size());
	for(size_t i=0; i<Go1.size();i++)
		{
		uG1(0,i)=Go1[i].x;
		uG1(1,i)=Go1[i].y;
		}
	for(size_t i=0; i<Go2.size();i++)
		{
		uG2(0,i)=Go2[i].x;
		uG2(1,i)=Go2[i].y;
		}
	// First Get the mean vectors
	ublas::vector<double> mG1;
	GetMeanByColumns<double>(uG1,mG1);
	ublas::vector<double> mG2;
	GetMeanByColumns<double>(uG2,mG2);
	//Get Meanof the mean vectors of each group
	// <G1>-<G2>
	ublas::vector<double> meandifG1G2=mG1-mG2;
	
	cout<<"Mean G1:"<<mG1<<endl;
	cout<<"Mean G2:"<<mG2<<endl;
	cout<<"Mean (G1-G2):"<<meandifG1G2<<endl;
	//substract mean  X=X-<X>
	//uG1=uG1-mG1;
	uG1=SubstractMeans(uG1, mG1);
	uG2=SubstractMeans(uG2, mG2);

	///////////////////////////
	matrix<double> uG1covar=prod(uG1,trans(uG1))/double(Go1.size());
	matrix<double> uG2covar=prod(uG2,trans(uG2))/double(Go2.size());

	matrix<double> pCovar=(double)(Go1.size()/(double)(Go1.size()+Go2.size()))*uG1covar + (double)(Go2.size()/(double)(Go1.size()+Go2.size()))*uG2covar;

	
	
	ublas::matrix<double> pCovarInv=pCovar;
	bool isinverted=InvertMatrix(pCovar, pCovarInv);
	ublas::vector<double> res=prod(pCovarInv, meandifG1G2);
	cout<<"Pooled Inverse Covar:"<<pCovarInv<<endl;
	cout<<"res=prod(pCovarinv, <G1,G2>):"<<res<<endl;
	//std::transform( res.begin(), res.end(), res.begin(),res.begin(),std::multiplies<double>() );// just pass the negated Mean
	double fres=std::accumulate(res.begin(), res.end(),0.0);
	cout<<"Covar:"<<pCovar<<endl;
	
	
	cout<<"Final result:"<<norm_2(res)<<endl;
	
	std::ifstream fin("c:/arm2arm/DATA/test_mahal.txt");
	float x, y;
	uG1.resize(2,1000);
	for(size_t i=0;i<1000;i++)
		{
		fin>>x>>y;
		uG1(0,i)=x;
		uG1(1,i)=y;
		}
	fin.close();
	Metrics::MahalDistance<double> Mahal(uG1);
	ublas::vector<double> vec(2);
	vec(0)=uG1(0,0);vec(1)=uG1(1,0);
	cout<<"Pos:"<<vec<<endl;
	cout<<"Dist: "<<Mahal.doDist(vec)<<endl;
	std::ofstream fout("c:/arm2arm/DATA/test_mahal1.txt");

	for(size_t i=0;i<1000;i++)
		{
		vec(0)=uG1(0,i);vec(1)=uG1(1,i);
		fout<<uG1(0,i)<<" "<<uG1(1,i)<<" "<<Mahal.doDist(vec)<<endl;
		}
	fout.close();
	
	return 0;

	}

