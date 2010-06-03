#ifndef _OAF_
#define _OAF_
#include "data_readers.h"

#define OAF32 32
#define OAF64 64
#define strOAF_32 "{32:fPOS3:fVEL3:fRHO:fHSML}"
#define strOAF_64 "{64:fPOS3:fVEL3:fRHO:fHSML:cPAD}"
#define strOAF_52 "{32:fPOS3:fVEL3:fMASS:fHSML:fRHO:fU:fZ2:fSFR}"

#define NCACHE 1024*1024

//int id, int snap, unsigned char type, float* pos, float *vel,
//			float mass, float hsml,float rho, float u, float *z, float sfr)
class COAFFile :public CGadget
	{

	public:
		COAFFile(std::string filename);
		COAFFile(void){};
		~COAFFile(void);
		bool Open(char io_flag='w');	
		bool Close();
		short GetFormat(){return m_Recformat;};
		bool WriteHeader(unsigned int Np,string strHead, double t1, double t2, double dt, int Ns);
		bool WriteHeader();
		bool WriteTime(char *pData, unsigned int timeSize);
		bool WriteOneParticle(char *buf, unsigned int Ns);
		unsigned int WriteDataStart(unsigned int Ns);
		bool WriteDataStop();

		private:

		ofstream m_ofile;
		char m_cacheOBUF[NCACHE];
		short m_Recformat;
		unsigned int m_datasize;
		unsigned int m_onerecsize;	
		struct strOAFhead{		
			string strHead;
			unsigned int Np;
			double Ts,Te,dT;
			unsigned int Ns;
			char *cPAD;
			}m_head;

	};

#endif
