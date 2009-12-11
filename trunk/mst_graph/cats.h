#ifndef MY_CATS
#define MY_CATS
#include <vector>
#include <iostream>
#include <fstream>
class TStr{
public:
	std::vector<uint>  id;
	uint m_Np;
	bool operator<(const TStr& that) const
		{
		if(m_Np < that.m_Np)
			return true;      
		return false;
		}
	};

class CCatalog{

public:
	std::vector<TStr> cats;//Catalogue

	void sort(){
		std::sort(cats.begin(),cats.end());
		std::reverse(cats.begin(),cats.end());
		}
	void SaveCats(std::string fname, uint nmax=10){
		std::ofstream of(fname.c_str(),std::ios::out | std::ios::binary);
		uint Nsave=std::min<int>(nmax, cats.size());
		of.write((char*)&Nsave, sizeof(Nsave));
		for(uint i=0;i<Nsave;i++)
			{
			of.write((char *)&cats[i].m_Np,sizeof(uint));
			of.write((char *)&cats[i].id[0],cats[i].m_Np*sizeof(uint));
			}
		of.close();
		}
	};

#endif