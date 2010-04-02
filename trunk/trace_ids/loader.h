#ifndef _LOADER_
#define _LOADER_

#include "data_readers.h"
#include "particleid.h"
#include <set>
template <class T>
class generator
	{
	T val_;
	T step_;
	public:
		generator(const T &val, const T step) : val_ (val),step_(step) { }
		T operator()() {
			return val_ += step_;
			}
	};
template<class T> T op_sum (T i, T j) { return i+j; }
template<class T> 
class strMover
	{
	T *pP_;
	T *pMove_;
	unsigned int stride_;
	public:
		strMover(T *pP,T *pMove, unsigned int stride):pP_(pP), pMove_(pMove),stride_(stride)
			{}
		void operator() ( const unsigned int ip )
			{
			std::transform(&pP_[ip*stride_], &pP_[ip*stride_]+stride_, 
				pMove_, &pP_[ip*stride_], op_sum<T>);
			}
	};

class CLoader
	{
	public:
		CLoader(std::string fname, std::set<int>  &IDfromFile, unsigned int ptype=0):m_fname(fname), m_ptype(ptype)
			{
				CGadget *pG=new CGadget(fname, false);
				m_nelem = pG->read_block<int>(pID,
                "ID  ",
                m_ptype);
				//m_nelem=10;/just for debugging
				if(m_nelem<1)
					{
					std::cout<<"Cannot open file, exiting"<<std::endl;
					exit(1);
					}
				if(IDfromFile.size()>0)
					{
					for(unsigned int i=0;i<m_nelem;i++)
						if(IDfromFile.find(pID[i])!=IDfromFile.end())
							data.insert(particleID(pID[i], i));
					}else
				for(unsigned int i=0;i<m_nelem;i++)
					data.insert(particleID(pID[i], i));
				///////////////////////////////////////////
				get_com_bypot(pG);
				unsigned int p_nelem = pG->read_blockv3(pPOS,
                "POS ",
                m_ptype);
				MoveToCOM();
				delete pG;
			}
		void get_com_bypot(CGadget *pG)
			{
				float *pP;
				unsigned int all_nelem = pG->read_block<float>(pP,
                "POT ",
                6);
				unsigned int comID=std::distance(pP, std::min_element(pP, pP+all_nelem));
				delete [] pP;
				unsigned int p_nelem = pG->read_blockv3(pPOS,
                "POS ",
                6);
				memcpy(&m_COM[0],&pP[comID*3],3*sizeof(float));
				delete [] pP;
				
			}
		
		void MoveToCOM(){
			unsigned int stride=3;
			for(unsigned int i=0;i<m_nelem;i++)
				{
					for(unsigned int j=0;j<stride;j++)
						pPOS[i*stride+j]+=m_COM[j];
				}
			generator<unsigned int> gen(0, 1);
			vector<unsigned int> idx(m_nelem);
			std::generate(idx.begin(), idx.end(),gen);
			std::for_each(idx.begin(), idx.end(), strMover<float>(&pPOS[0],&m_COM[0],stride));
			}
		void PrintStat(){
			if(m_nelem>0)
				{
				cout<<"Minval:"<<(*std::min_element(pID, pID+m_nelem))<<endl;
				cout<<"Maxval:"<<(*std::max_element(pID, pID+m_nelem))<<endl;
				}
			};
		~CLoader()
			{
			if(pID!=NULL)
				delete[] pID;
			if(pPOS!=NULL)
				delete [] pPOS;
			std::cout<<"exiting..."<<m_fname<<std::endl;
			}
		particlesID_set data;//keeping relation for ID, IDf
		std::string m_fname;
		unsigned int m_nelem;
		unsigned int m_ptype;
		int *pID;
		float *pPOS;
		float m_COM[3];
	};


#endif