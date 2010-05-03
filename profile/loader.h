#ifndef _LOADER_
#define _LOADER_

#include "data_readers.h"
#include "particleid.h"
#include "utils.h"
#include <functional>//negate
#include <algorithm>//transform

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
		CLoader(std::string fname, unsigned int ptype=0):m_fname(fname),pID(NULL), pType(NULL),m_ptype(ptype)
			{
				CGadget *pG=new CGadget(fname, false);
pG->m_verbose=0;
				//m_nelem = pG->read_block<int>(pID,
                //"ID  ",
                //m_ptype);
		//		//m_nelem=10;/just for debugging
		//		if(m_nelem<1)
		//			{
		//			std::cout<<"Cannot open file, exiting"<<std::endl;
		///			exit(1);
		//			}
				/*for(unsigned int i=0;i<m_nelem;i++)
					data.insert(particleID(pID[i], i));*/
				///////////////////////////////////////////
				get_com_bypot(pG);
                                
				unsigned int p_nelem = pG->read_blockv3(pPOS,
                "POS ",
                m_ptype);
                                 m_nelem=p_nelem;
				if(m_ptype==6)
					{
					pType= new int[m_nelem];
					for(int it=0,ip=0;it<6;it++)
						for(int i=0;i<pG->myhead.npart[it];i++)
							{
							  pType[ip++]=it;
							}
					int np=pG->myhead.npart[0]+pG->myhead.npart[4];
					int npout=pG->read_full_block<float>(pSFR, "SFR ", np);
					if(np!=npout){cout<<"error reading SFR"<<endl;exit(0);} 
					}
				memcpy(&m_npart[0], &pG->myhead.npart[0], 6*sizeof(unsigned int));
				MoveToCOM();
				delete pG;
			}
		float R(unsigned int i){
			float val=pPOS[i*3+0]*pPOS[i*3+0]+
				pPOS[i*3+1]*pPOS[i*3+1]+
				pPOS[i*3+2]*pPOS[i*3+2];
			return sqrt(val);
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
			
			std::transform ( m_COM, m_COM+3, m_COM, std::negate<float>() );
		
			generator<unsigned int> gen(0, 1);
			vector<unsigned int> idx(m_nelem);
			std::generate(idx.begin(), idx.end(),gen);
			std::for_each(idx.begin(), idx.end(), strMover<float>(&pPOS[0],&m_COM[0],stride));
			std::transform ( m_COM, m_COM+3, m_COM, std::negate<float>());
			}
		void PrintStat(){
			if(m_nelem>0)
				{
				cout<<"\n============ ID ============= "<<endl;
				/*cout<<"Minval:"<<(*std::min_element(pID, pID+m_nelem))<<endl;
				cout<<"Maxval:"<<(*std::max_element(pID, pID+m_nelem))<<endl;
*/
				cout<<"\n============ SFR ============= "<<endl;
				cout<<"Minval:"<<(*std::min_element(pSFR, pSFR+m_npart[0]+m_npart[4]))<<endl;
				cout<<"Maxval:"<<(*std::max_element(pSFR, pSFR+m_npart[0]+m_npart[4]))<<endl;
				}
			};
		void ReadAP(std::string file)
			{
			std::ifstream infile(file.c_str());

			int val;
			string f1, f2;
			float A, P, AP;
			if(is_file_exist(file))
				{
				//infile>>f1;
				//infile>>f2;

				while(infile>>val>>A>>P>>AP)
					{
//					m_IDlist.insert(std::make_pair(val,TApData(A,P)));
					}
				}
			else
				std::cerr<<"WARNING::Cannot find ID file: "<<file<<std::endl;
			};
		~CLoader()
			{
			if(pID!=NULL)
				delete[] pID;
			if(pType!=NULL)
				delete[] pType;
			if(pSFR!=NULL)
				delete[] pSFR;
			if(pPOS!=NULL)
				delete [] pPOS;
			std::cout<<"exiting..."<<m_fname<<std::endl;
			}
		particlesID_set data;//keeping relation for ID, IDf
		std::string m_fname;
		unsigned int m_nelem;
		unsigned int m_npart[6];
		unsigned int m_ptype;
		int *pID;
		int *pType;
		float *pSFR;
		float *pPOS;
		float m_COM[3];
	};


#endif




