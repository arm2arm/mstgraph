#ifndef _MYSIGMA_
#define _MYSIGMA_
////////////// This class does velosity disperson profile.
//usage: 	
// CSigma<double> GetSigma(data.x, data.y, data.z, data.vx, data.vy, data.vz);
//
#include <vector>      // for vector
#include <algorithm>   // for adjacent_find, find
#include <functional>  // for bind1st, equal_to
#include <iostream>    // for cout, endl
#include <fstream>    // for file IO
#include <stdlib.h>
#include <cmath>    // for sqrt
#include <cstring>    // for string
#include <string>    // for string
#include <bitset> //bitsets

#include <valarray> //vallaray
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
//#include "utils.h"
using std::vector;
using std::valarray;
using std::cout;
using std::cerr;
using std::endl;
using std::string;


/////////////// UTILS for some mnipulation
struct TWhereIs{
	vector<int> idx;
	vector<bool> hash;
	TWhereIs(vector<bool> bv){
		hash=bv;	
		};
	TWhereIs(){};
	vector<int> get(){return idx;};
	void update(){
		idx.clear();
		for(size_t i=0;i<hash.size();i++)
			if(hash[i])
				idx.push_back(i);
		};
	void reset(){hash.clear();idx.clear();};
	void push(vector<bool> bv){
		if(hash.size()==0){
			hash=bv;			
			}else{
				if(bv.size() != hash.size())
					cerr<<"# Error in TWhereIs the vectors are not equal."<<endl;
				typedef vector<bool> BV;
				BV::iterator ibhash=hash.begin(), iehash=hash.end();
				for(int i=0;iehash!=ibhash;ibhash++){
					(*ibhash) = (*ibhash) && bv[i++];
					}
			}
		update();
		};
	inline size_t size(){
		return (int)std::count(hash.begin(), hash.end(), true);
		};
	};

template <class Tg>
valarray<Tg> get(vector<int> &ib,vector<Tg> &v)
	{
	valarray<Tg> tmp(ib.size());
	for(size_t i=0;i<ib.size();i++)
		{
		tmp[i]= v[ib[i]];
		}
	return tmp;
	};


void mma(vector<double> &data, std::string msg="Stat ")
	{

	std::cout<<msg+string(" min=")<<(*std::min_element(data.begin(), data.end())) <<endl;
	std::cout<<msg+string(" max=") << *std::max_element(data.begin(), data.end())<<endl;
	};



////////////////////////////////////////////
template <class T>
class CSigma{
public:
	vector< vector<double> > sigma; 
	vector<double> rr, rrslw;
	string m_fname;
	void save(vector<double> &rr_, vector< vector<double> > &sig, string fname)
		{
		std::ofstream stream(fname.c_str());
		if(stream.is_open())
			{
			stream<<"# Rad \t";
			for(size_t ity=0;ity<types.size();ity++)
				stream<<"sigma_T"<<types[ity]<<"\t";
			stream<<endl;
			for(size_t i=0;i<rr_.size();i++)
				{
				stream<<rr_[i]<<"\t";
				for(size_t ity=0;ity<types.size();ity++)
					stream<<sig[ity][i]<<"\t";
				stream<<endl;
				}
			}
		}
	struct CData
		{
		inline size_t size(){return x.size();};
		inline T R(size_t i)
			{
			if(i<x.size())
				return sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
			return 0;
			}
		inline void reserve(size_t np){
			m_np=np;
			x.resize(np);
			y.resize(np);
			z.resize(np);
			vx.resize(np);
			vy.resize(np);
			vz.resize(np);
			type.resize(np);
			};
		void insert(int type_,float *pX, float *pV)
			{
			//	if(abs(pX[2])<3.0 && abs(pV[2])<400.0 )
					{
				x.push_back(pX[0]);
				y.push_back(pX[1]);
				z.push_back(pX[2]);
				vx.push_back(pV[0]);
				vy.push_back(pV[1]);
				vz.push_back(pV[2]);
				type.push_back(type_);
				}
			}
		void insert(int type_,T x_, T y_,T  z_, T vx_, T vy_, T vz_)
			{
				x.push_back(x_);
				y.push_back(y_);
				z.push_back(z_);
				vx.push_back(vx_);
				vy.push_back(vy_);
				vz.push_back(vz_);
				type.push_back(type_);
			}
		valarray<double> GetDistXY()
			{
			valarray<double> dist(this->size());
			for(size_t i=0;i<this->size();i++)
				{
				dist[i]=sqrt(x[i]*x[i]+y[i]*y[i]);
				}
			return  dist;
			}
		
		size_t m_np;
		vector<T> x,y,z,vx, vy, vz;
		vector<int> type;
		void SubstractMean(vector<T> &x)
			{
			MeanValue mv = for_each (x.begin(), x.end(),  // range
				MeanValue());              // operation
			std::transform(x.begin(), x.end(), x.begin(), std::bind2nd(std::plus<T>(), -mv.value()));
			};
		void PutInCom(void)
		{
			
		SubstractMean(x);SubstractMean(y);SubstractMean(z);
		SubstractMean(vx);SubstractMean(vy);SubstractMean(vz);

		};
		}data;
		enum eTYPE{T0,T1,T2,T3,T4,T5,numtypes};
		std::bitset<numtypes> userTypes;
		vector<int> types;
		void GetTypes(std::string &strType){
			string  names[]={"T0","T1","T2","T3","T4","T5"};
			std::string::iterator its=strType.begin();
			while( its!=strType.end() )
				{
				char ch=*its;
				int val=atoi(&ch);
				if(val > -1 && val<numtypes+1)
					{
					userTypes.set(val);
					types.push_back(val);
					}
				cout<<*its<<endl;
				its++;
				}
			};

	CSigma(std::string strType,int *pType, float *pX, float *pV,size_t np)
		{
		m_fname="sigma.txt";
		GetTypes(strType);
		if (userTypes.any())
			{
			for(size_t i=0;i<np;i++)
				{
				if (userTypes[(eTYPE)pType[i]])
					{
					data.insert(pType[i],&pX[i*3],&pV[i*3]);
					}
				}
			
			cout<<"# Prepare for sigma over the "<<data.size()<<" particles"<<endl;

			}
		};
	~CSigma(){

		save(rr,sigma, m_fname);
		save(rrslw, sigmaSlitX, m_fname+".slitX");
		save(rrslw, sigmaSlitY, m_fname+".slitY");
		}

	vector<double> Jv;
	vector<double> GetJv(void)
		{
		Jv.resize(3);
		size_t np=0;
		for(size_t i =0; i<data.vz.size();i++)
			{
			if(data.type[i]==4 && data.R(i)<40.0)
				{
				Jv[0] += (data.y[i]*data.vz[i]-data.z[i]*data.vy[i]);
				Jv[1] += (data.z[i]*data.vx[i]-data.x[i]*data.vz[i]);
				Jv[2] += (data.x[i]*data.vy[i]-data.y[i]*data.vx[i]);
				np++;
				}
			}
		Jv[0]/=static_cast<double>(np);
		Jv[1]/=static_cast<double>(np);
		Jv[2]/=static_cast<double>(np);
		return Jv;
		}
	void GetSigma(size_t Nbins=150, double Rc=4.0){
		

		int typeCount=userTypes.count();
		if(!typeCount)return;
		
		cout<<"# Sigma in the X and Y slits...";
		cout.flush();
		GetSigmaInSlit(Nbins, Rc);
		cout<<"..done"<<endl;
		
		
		cout<<"# Sigma within circular rings ..";
		cout.flush();
		valarray<double> dist(0.0,data.size());
		for(size_t i=0;i<data.size();i++)
			{
			dist[i]=sqrt(data.x[i]*data.x[i]+data.y[i]*data.y[i]);
			}
		
		double dr=0.1;
		Nbins=(size_t)(Rc/dr);
		double r;
		sigma.resize(6);
		for(size_t i=0;i<6;i++)
		sigma[i].resize(Nbins);
		rr.resize(Nbins);
		valarray<double> vz(&data.vz[0], data.vz.size());
		//valarray<double> z(&data.z[0], data.z.size());
		
		valarray<int>    type(&data.type[0], data.type.size());

		double mvel=vz.sum()/(double)vz.size();
		for(size_t i=0l;i<Nbins;i++)
			{
			r=dr*i;
			rr[i]=r;
			//cout<<i<<") "<<r<<" ";
			for(size_t itype=0;itype<types.size();itype++)
				{
				valarray<bool> ids = (dist < r+dr) && (dist > r) && (type==types[itype]);
				if(ids.max())
					{
					valarray<double> d=vz[ids];
					d-=mvel;
					d=pow(d,2.0);
					sigma[itype][i]= sqrt(d.sum()/(double)d.size());
					}
				//cout<<sigma[itype][i]<<"\t";
				}
			//cout<<endl;
			}
		cout<<"..done"<<endl;
		for(size_t i=0;i<6;i++)
			smooth(sigma[i]);
		
		};
/////////////////////////////// 
	vector<vector<double> > sigmaSlitX, sigmaSlitY;
	void GetSigmaInSlit(size_t Nbins=150, double Rc=4.0, double slw=0.5, double dr=0.1){
		CData slitdataX, slitdataY;
		if(false)
			{std::ofstream stream("test.txt"),streamy("test2.txt");
		if(stream.is_open());
			}
		for(size_t i=0;i<data.size();i++)
			{
			if(abs(data.y[i])<slw && abs(data.x[i])<Rc)
				{
				slitdataX.insert(data.type[i], 
				data.x[i], data.y[i],  data.z[i],  data.vx[i],  data.vy[i],  data.vz[i] );
		if(false)		stream<<data.x[i]<<" "<<data.y[i]<<endl;
				}
			if(abs(data.x[i])<slw && abs(data.y[i])<Rc)
				{slitdataY.insert(data.type[i], 
				data.x[i], data.y[i],  data.z[i],  data.vx[i],  data.vy[i],  data.vz[i] );
			if(false)streamy<<data.x[i]<<" "<<data.y[i]<<endl;
				}
			}
		if(false){stream.close();streamy.close()};
		

		dr=(2.0*Rc)/static_cast<double>(Nbins);
		cout<<"Nbins="<<Nbins<<endl;
		double r;
		sigmaSlitX.resize(6);
		sigmaSlitY.resize(6);
		for(size_t i=0;i<6;i++)
			{
			sigmaSlitX[i].resize(Nbins);
			sigmaSlitY[i].resize(Nbins);
			}
		rrslw.resize(Nbins);
		
		valarray<double> slitXvz(&slitdataX.vz[0], slitdataX.vz.size());
		valarray<int>    typeX(&slitdataX.type[0], slitdataX.type.size());
		
		valarray<double> slitYvz(&slitdataY.vz[0], slitdataY.vz.size());
		valarray<int>    typeY(&slitdataY.type[0], slitdataY.type.size());

		double mvelX=slitXvz.sum()/static_cast<double>(slitXvz.size());
		double mvelY=slitYvz.sum()/static_cast<double>(slitYvz.size());

		valarray<double> slwX(&slitdataX.x[0], slitdataX.x.size());		
		valarray<double> slwY(&slitdataY.y[0], slitdataY.y.size());		

		for(size_t i=0l;i<Nbins;i++)
			{
			r=dr*i-Rc;
			rrslw[i]=r;
			if(false)cout<<i<<" "<<r<<endl;
			for(size_t itype=0;itype<types.size();itype++)
				{
				valarray<bool> idsX = (slwX < r+dr) && (slwX > r) && (typeX==types[itype]);
				valarray<bool> idsY = (slwY < r+dr) && (slwY > r) && (typeY==types[itype]);	
				if(idsX.max())
					sigmaSlitX[itype][i]=check_and_get(idsX, mvelX,slitXvz);
				if(idsY.max())
					sigmaSlitY[itype][i]=check_and_get(idsY, mvelY,slitYvz);
				}
			}
		for(size_t i=0;i<6;i++)
			{
			smooth(sigmaSlitX[i]);
			smooth(sigmaSlitY[i]);
			}
		//save(rrslw, sigmaSlitX, m_fname+".slitX");
		//save(rrslw, sigmaSlitY, m_fname+".slitY");
		
		};
///////////////////////////////
	double check_and_get(valarray<bool> &ids_, double mvel/*meanvalue of the velociti in the whole region*/, 
		valarray<double> &vz_)
		{
		double sig=0.0;
		if(ids_.max())
			{
			valarray<double> d=vz_[ids_];
			d-=mvel;
			d=pow(d,2.0);
			sig = sqrt(d.sum()/(double)d.size());
			}
		return sig;
		}

	};

#endif
