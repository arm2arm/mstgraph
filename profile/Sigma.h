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

#include <cmath>    // for sqrt
#include <cstring>    // for string
#include <string>    // for string

#include <valarray> //vallaray
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
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
	vector<int>::iterator it=ib.begin();
	for(int i=0;it!=ib.end();it++)
		{
		tmp[i]= v[(*it)];
		i++;
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
	vector<double> sigma; 
	vector<double> rr;
	string m_fname;
	void save()
		{
		std::ofstream stream(m_fname.c_str());
		if(stream.is_open())
			for(size_t i=0;i<rr.size();i++)
				{
				stream<<rr[i]<<" "<<sigma[i]<<endl;
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
			if(i<m_np)
				{
				x[i]=pX[0];
				y[i]=pX[1];
				z[i]=pX[2];
				vx[i]=pV[0];
				vy[i]=pV[1];
				vz[i]=pV[2];
				type[i]=type_;
				}
			}
		size_t m_np;
		vector<T> x,y,z,vx, vy, vz;
		vector<int> type;
		}data;

	CSigma(int *pType, float *pX, float *pV,size_t np)
		{
		m_fname="sigma.txt";
		data.reserve(np);	
		for(size_t i=0;i<np;i++)
			{
			data.insert(i,pType[i],&pX[i*3],&pV[i*3]);
			}
		GetSigma(&data, 6, sigma, rr, 10);
		};
	~CSigma(){

		save();
		}
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


	void GetSigma(CData *pL,int type, vector<double> &sigma, vector<double> &rr, double Rc=5.0)
		{

		size_t Nsigbin=150, Nbins=100;
		double dr=Rc/(double)Nsigbin;
		double slw=Rc;
		vector<int> ids,ib;
		sigma.resize(Nsigbin, (double)0.0);
		rr=sigma;

		vector<double> dist;
		for(size_t i=0;i<pL->size(); i++)
			dist.push_back(pL->R(i));

		TWhereIs whereib;
		//check (x> -rc && x<rc && abs(y)< slw)
		vector<bool> ans=make_bool_vec<double, double>(pL->x, -Rc, std::greater_equal<double>()); 
		whereib.push(ans);
		ans=make_bool_vec<double, double>(pL->x, Rc, std::less<double>()); 
		whereib.push(ans);
		ans=make_bool_vec<double, double>(pL->y, -Rc, std::greater_equal<double>()); 
		whereib.push(ans);
		ans=make_bool_vec<double, double>(pL->y, Rc, std::less<double>()); 
		whereib.push(ans);

		ans=make_bool_vec<int, int>(pL->type, 4, std::less<int>()); 
		whereib.push(ans);

		ib=whereib.get();
		mma(dist, "Dist: ");
		TWhereIs whereis;	
		for(size_t i=0;i<Nsigbin;i++)
			{
			double r=i*dr;
			whereis.reset();
			whereis.push(
				make_bool_vec<double, double>(dist, r, std::greater_equal<double>())
				);
			whereis.push(
				make_bool_vec<double, double>(dist, r+dr, std::less<double>())
				);
			cout<<r<<endl;
			ids=whereis.get();
			if(ids.size() > 0)
				{
				valarray<double> vz=get(ib,pL->vz);
				double meanvelz=vz.sum()/(double)vz.size();
				vz=get(ids, pL->vz);
				vz-=meanvelz;
				vz=std::pow(vz,(double)2.0);
				sigma[i]=sqrt(vz.sum()/(double)ids.size());
				rr[i]=r;
				}
			}

		}
	};

#endif
