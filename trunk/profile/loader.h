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
			std::transform(&(pP_[ip*stride_]), &(pP_[ip*stride_ +stride_]), 
				pMove_, &(pP_[ip*stride_]), op_sum<T>);
			}
	};

class CLoader
	{
	public:
		CLoader(std::string fname, unsigned int ptype=0, bool velflag=false):
		  m_fname(fname),
			  pSFR(NULL),
			  pID(NULL), 
			  pType(NULL),pVEL(NULL),pMASS(NULL),pU(NULL),pRHO(NULL),pZ(NULL),pHSML(NULL),m_ptype(ptype), m_verbose(0)
			  {
			  CGadget *pG=new CGadget(fname, false);
			  pG->m_verbose=0;
			  if(m_ptype==5&&pG->myhead.npart[5]<100)
				  {
				  //     pG->m_verbose=2;
				  int np=pG->myhead.npart[5];
				  float *tdata=NULL;//=new float[np];
				  unsigned int nelem = pG->read_full_block<float>(tdata, "BHMA", np);
				  if(nelem>0)
					  m_BHDATA[0]=tdata[0];
				  nelem = pG->read_full_block<float>(tdata, "BHMD", np);
				  if(nelem>0)
					  m_BHDATA[1]=tdata[0];
				  if(tdata!=NULL)delete [] tdata;

				  } 

			  unsigned int p_nelem = pG->read_blockv3(pPOS,"POS ", m_ptype);
			  m_nelem=p_nelem;
			  if(velflag)
				  pG->read_blockv3(pVEL,"VEL ", m_ptype);
			  if(m_ptype==6)
				  {
				  pType= new int[m_nelem];
				  for(int it=0,ip=0;it<6;it++)
					  for(int i=0;i<pG->myhead.npart[it];i++)
						  {
						  pType[ip++]=it;
						  }
					  int np=pG->myhead.npart[0];
					  int npout=pG->read_full_block<float>(pSFR, "SFR ", np);
					  if(np!=npout){cout<<"error reading SFR"<<endl;exit(0);}
					  npout=pG->read_full_block<float>(pU, "U   ", np);
					  npout=pG->read_full_block<float>(pRHO, "RHO ", np);
					  npout=pG->read_full_block<float>(pHSML, "HSML", np);
					  np=pG->myhead.npart[0]+pG->myhead.npart[4];
					  npout=pG->read_full_block<float>(pZ, "Z   ", np);
					  npout=pG->read_full_block<float>(pMASS, "MASS", np);
					  
				  }
			  memcpy(&m_npart[0], &pG->myhead.npart[0], 6*sizeof(unsigned int));
			  
			  delete pG;
			  }
		  inline 		float R(unsigned int i){
			  float val=pPOS[i*3+0]*pPOS[i*3+0]+
				  pPOS[i*3+1]*pPOS[i*3+1]+
				  pPOS[i*3+2]*pPOS[i*3+2];
			  return sqrt(val);
			  }
		  inline 		float Get3DVel(unsigned int i){
			  float val=pVEL[i*3+0]*pVEL[i*3+0]+
				  pVEL[i*3+1]*pVEL[i*3+1]+
				  pVEL[i*3+2]*pVEL[i*3+2];
			  return sqrt(val);
			  }
		  inline 		float GetBHMass(int i=0){return m_BHDATA[i*2];}
		  inline 		float GetBHDOTMass(int i=0){return m_BHDATA[i*2+1];}
		  void get_com_bypot(CGadget *pG)
			  {
			  float *pP;
			  unsigned int all_nelem = pG->read_block<float>(pP,"POT ", 6);
			  unsigned int comID=std::distance(pP, std::min_element(pP, pP+all_nelem));
			  delete [] pP;
			  unsigned int p_nelem = pG->read_blockv3(pPOS,"POS ", 6);
			  memcpy(&m_COM[0],&pPOS[comID*3],3*sizeof(float));
			  delete [] pPOS;

			  }

		  template<class T>
		  void MoveToCOM(T *com=NULL,unsigned int ncom=3 ){
			  const unsigned int stride=3;
			  if(com!=NULL)
				  for(size_t i=0;i<ncom;i++)
					  m_COM[i]=static_cast<float>(com[i]);
			  for(size_t ip=0;ip<m_nelem;ip++)
				  {
				  pPOS[ip*stride]   -= m_COM[0];
				  pPOS[ip*stride+1] -= m_COM[1];
				  pPOS[ip*stride+2] -= m_COM[2]; 
				  }
			  if(ncom ==6 && pVEL!=NULL)
				  for(size_t ip=0;ip<m_nelem;ip++)
					  {
					  pVEL[ip*stride]   -= m_COM[3];
					  pVEL[ip*stride+1] -= m_COM[4];
					  pVEL[ip*stride+2] -= m_COM[5]; 
					  }

				  /* std::transform ( m_COM, m_COM+3, m_COM, std::negate<float>() );		  
				  generator<unsigned int> gen(0, 0);
				  vector<unsigned int> idx(m_nelem);
				  std::generate(idx.begin(), idx.end(),gen);
				  std::for_each(idx.begin(), idx.end(), strMover<float>(&pPOS[0],&m_COM[0],stride));
				  std::transform ( m_COM, m_COM+3, m_COM, std::negate<float>());*/
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
		  bool ReadID(){
			  CGadget *pG=new CGadget(m_fname, false);
			  pG->m_verbose=0;
			  unsigned int np=pG->read_block<int>(pID, "ID  ", m_ptype);
			  if(np!=size())return false;
			  return true;
			  delete pG;
			  }
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
			  if(pVEL!=NULL)
				  delete [] pVEL;
  if(pU!=NULL)
				  delete [] pU;
			  if(pZ!=NULL)
				  delete [] pZ;
			  if(pHSML!=NULL)
				  delete [] pHSML;
			  if(pRHO!=NULL)
				  delete [] pRHO;
			  if(pMASS!=NULL)
				  delete [] pMASS;
						
			  if(m_verbose==2)std::cout<<"exiting..."<<m_fname<<std::endl;
			  }
		  inline size_t size(){return m_nelem;}
		  particlesID_set data;//keeping relation for ID, IDf
		  std::string m_fname;
		  size_t m_nelem;
		  unsigned int m_npart[6];
		  unsigned int m_ptype;
		  int *pID;
		  int *pType;
		  float *pSFR;
		  float *pU;
		  float *pRHO;
		  float *pZ;
		  float *pMASS;
		  float *pHSML;
		  float *pPOS;
		  float *pVEL;
		  float m_COM[6];
		  float m_BHDATA[2];//Mbh, DOTMbh
		  int m_verbose;
	};


#endif




