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
	vector<double> rr;
	string m_fname;
	void save()
		{
		std::ofstream stream(m_fname.c_str());
		if(stream.is_open())
			{
			stream<<"# Rad \t";
			for(size_t ity=0;ity<types.size();ity++)
				stream<<"sigma_T"<<types[ity]<<"\t";
			stream<<endl;
			for(size_t i=0;i<rr.size();i++)
				{
				stream<<rr[i]<<"\t";
				for(size_t ity=0;ity<types.size();ity++)
					stream<<sigma[ity][i]<<"\t";
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
		void insert(size_t i,int type_,float *pX, float *pV)
			{
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
					data.insert(i,pType[i],&pX[i*3],&pV[i*3]);
					}
				}
			
			cout<<"# Sigma over the "<<data.size()<<" particles"<<endl;

			}
		};
	~CSigma(){

		save();
		}
	void GetSigma(size_t Nbins=150, double Rc=4.0){

		int typeCount=userTypes.count();
		if(!typeCount)return;

		valarray<double> dist(0.0,data.size());
		for(size_t i=0;i<data.size();i++)
			{
			dist[i]=sqrt(data.x[i]*data.x[i]+data.y[i]*data.y[i]);
			}
		
		double dr=Rc/(double)Nbins, r;
		sigma.resize(6);
		for(size_t i=0;i<6;i++)
		sigma[i].resize(Nbins);
		rr.resize(Nbins);
		valarray<double> vz(&data.vz[0], data.vz.size());
		valarray<int>    type(&data.type[0], data.type.size());

		double mvel=vz.sum()/(double)vz.size();
		for(size_t i=0l;i<Nbins;i++)
			{
			r=dr*i;
			cout<<i<<") "<<r<<" ";
			for(size_t itype=0;itype<types.size();itype++)
				{
				valarray<bool> ids = (dist < r+dr) && (dist > r) && (type==types[itype]);
				if(ids.max())
					{
					valarray<double> d=vz[ids];
					d-=mvel;
					d=pow(d,2.0);
					sigma[itype][i]= sqrt(d.sum()/(double)d.size());
					rr[i]=r;
					}
				cout<<sigma[itype][i]<<"\t";
				}
			cout<<endl;
			}
		};
	template <class Tvec, class Tconst, typename TOpbin>
	vector<bool> make_bool_vec(vector<Tvec> &vec, Tconst value, TOpbin op, int *np=NULL)
		{
		vector<bool> myanswer(vec.size(), false);	
		typename vector<Tvec>::iterator ie=vec.end();
		typename vector<Tvec>::iterator ib=vec.begin();
		size_t i=0;
		for(i=0;i<vec.size();i++)
			{
			if( op((*ib), value))
				{
				myanswer[i].flip();
				}
			ib++;
			}
		int nc=(int)std::count(myanswer.begin(), myanswer.end(), true);
		if(np!=NULL)(*np)=i;
		return myanswer;
		};


	void GetSigma(CData *pL, vector<double> &sigma, vector<double> &rr, double Rc=5.0)
		{
		size_t Nsigbin=250, Nbins=200;
		double dr=Rc/(double)Nsigbin;
		double slw=Rc;
		vector<int> ids,ib;
		sigma.resize(Nsigbin, (double)0.0);
		rr=sigma;

		vector<double> dist;
		for(size_t i=0;i<pL->size(); i++)
			dist.push_back(sqrt(pL->x[i]*pL->x[i]+pL->y[i]*pL->y[i]));
                cout<<dist.size()<<endl;
		TWhereIs whereib;
		//check (x> -rc && x<rc && abs(y)< slw)
		vector<bool> ans=make_bool_vec<double, double>(dist, Rc, std::less_equal<double>()); 
		whereib.push(ans);
		/*vector<bool> ans=make_bool_vec<double, double>(pL->x, -Rc, std::greater_equal<double>()); 
		whereib.push(ans);
		ans=make_bool_vec<double, double>(pL->x, Rc, std::less<double>()); 
		whereib.push(ans);
		ans=make_bool_vec<double, double>(pL->y, -Rc, std::greater_equal<double>()); 
		whereib.push(ans);
		ans=make_bool_vec<double, double>(pL->y, Rc, std::less<double>()); 
		whereib.push(ans);*/
		ib=whereib.get();
		mma(dist, "Dist: ");
		valarray<double> vzb=get(ib,pL->vz);
		double meanvelz=vzb.sum()/(double)vzb.size();

		for(size_t i=0;i<Nsigbin;i++)
			{
			double r=i*dr;
			TWhereIs whereis;	
			whereis.push(
				make_bool_vec<double, double>(dist, r, std::greater_equal<double>())
				);
			whereis.push(
				make_bool_vec<double, double>(dist, r+dr, std::less<double>())
				);
			ids=whereis.get();
			if(ids.size() > 0 )
				{
				valarray<double> vzslice=get(ids,pL->vz);
				vzslice=std::pow(vzslice-meanvelz,(double)2.0);
				sigma[i]=sqrt(vzslice.sum()/(double)ids.size());
				rr[i]=r;
				}
			}

		}
	};

#endif
