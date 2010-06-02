#include "oaffile.h"



COAFFile::COAFFile(std::string filename)
	{
	m_filename=filename+string(".oaf.0");
	cout<<"# COAFFile::COAFFile will  write to "<<m_filename<<endl;

	}

COAFFile::~COAFFile(void)
	{
	}
bool COAFFile::Open(char io_flag)
	{
	bool ret_flag=true;
	if(io_flag=='r')
		{
		std::cerr<<"0 not implemented"<<endl;ret_flag=false;
		}else
		{
		swp_flag=false;
		m_ofile.open(m_filename.c_str(),std::ios::binary);
		m_ofile.rdbuf()->pubsetbuf(m_cacheOBUF, NCACHE);
			}

		if(!m_ofile.good())
			{
			cerr<<"Fail to open file... "<<m_filename<<"\nExiting..."<<endl;exit(0);
			}else  cerr<<"File is ok... "<<endl;

		return ret_flag;
	}
bool COAFFile::WriteHeader()
	{
	return WriteHeader(10,string(strOAF_32), 1.0, 10.0, 0.2, int(9.0/0.2));

	}
bool COAFFile::WriteHeader(unsigned int Np,string strHead, double t1, double t2, double dt, int Ns)
	{

	string tstr(strHead.begin()+1, strHead.begin()+ strHead.find_first_of(":"));;

	m_Recformat = (short)atol(tstr.c_str());

	m_head.strHead=strHead;
	m_head.Ns=Ns;
	m_head.Np=Np;
	m_head.Ts=t1;
	m_head.Te=t2;
	m_head.dT=dt;
	//  unsigned int shead=sizeof(m_head);


	char *buf=(char *)malloc(1024);
	memset(buf, 'X', 1024);
	strcpy(buf, m_head.strHead.c_str());
	memcpy((void*)&buf[0+m_head.strHead.length()], (void*)&m_head.Np, 4);
	memcpy((void*)&buf[4+m_head.strHead.length()], (void*)&m_head.Ts, 8);
	memcpy((void*)&buf[12+m_head.strHead.length()], (void*)&m_head.Te, 8);
	memcpy((void*)&buf[20+m_head.strHead.length()], (void*)&m_head.dT, 8);
	memcpy((void*)&buf[28+m_head.strHead.length()], (void*)&m_head.Ns, 4);

	WriteOneBlock(m_ofile,"HEAD", 
		(char*)&buf[0], 1024);
	free(buf);
	if(!m_ofile.good())
		{
		cerr<<"Fail to write HEAD record name to file...\nExiting..."<<endl;
		}
	return true;

	}
bool COAFFile::WriteTime(char *pData, unsigned int timeSize)
	{
	WriteOneBlock(m_ofile,"TIME", 
		(char*)&pData[0], timeSize);	 		
	return true;
	}

bool COAFFile::WriteDataStop()
	{
	m_ofile.write((char*)&m_datasize,4);	
	return true;
	}

unsigned int  COAFFile::WriteDataStart(unsigned int Ns)
	{
	if(Ns!=m_head.Ns) {cout<<"strange error"<<endl;
	return false;}
	unsigned int blsize, idata;
	m_onerecsize=(m_Recformat*m_head.Ns+4);
	m_datasize=m_onerecsize*m_head.Np;
	/*write block name*/
	blsize=8;
	m_ofile.write((char*)&blsize,4);
	m_ofile.write("DATA",4);		
	idata=8+m_datasize;
	m_ofile.write((char*)&idata,4);
	m_ofile.write((char*)&blsize,4);

	///Write F77 starting block ///
	m_ofile.write((char*)&idata,4);
	if(!m_ofile.good())
		{
		cerr<<"Fail to write DATA record name to file...\nExiting..."<<endl;
		}
	return m_onerecsize;
	}

bool COAFFile::WriteOneParticle(char *pData, unsigned int indatasize)
	{
	//  static long long iglobal_data;
	if(indatasize!=(m_onerecsize)){cout<<"Wrong datasize"<<endl;exit(1);}
	m_ofile.write((char*)pData,indatasize);

	if(!m_ofile.good())
		{
		cerr<<"Fail to write ID record to file...\nExiting..."<<endl;
		cerr<<"Size was : "<<indatasize<<endl;
		exit(1);
		}
	return  true;
	}
bool COAFFile::Close()
	{ 			
	m_ofile.close();
	return true;
	}
