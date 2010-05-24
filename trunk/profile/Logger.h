#ifndef _LOGGER_
#define _LOGGER_
#include "loader.h"
#include "utils.h"
#include <iterator>
struct TData
	{
	std::vector<float> v;
	size_t size(){return v.size();};
	friend std::ostream & operator<<(std::ostream& , const TData& p );
	};
std::ostream & operator << (std::ostream &os, const TData &p)
	{
	for(size_t i=0;i<p.v.size();i++)
	os << p.v[i]<<"\t";
	return os;
	}

class CLogger :public CAsciiReader{
public:
	CLogger(std::string fout="dump.log", bool append=false,int nfields=1):m_write_on_save(true),m_append(append),CAsciiReader(fout, nfields,append ){
				    
		};
	~CLogger(){if(m_write_on_save)flush();};
	void flush(){save();};
	////////// Populate ////////
	
	bool insert(const int isnap, const dynvector &d)
		{
		bool retval=is_done(isnap);
		m_data[isnap]=d;
		flush();
		return retval;
		
		};
        bool is_done(unsigned int isnap)
	{return m_append&&(m_data.find(isnap)!=m_data.end());}
	//////////// IO ////////
	//		void load()//moved to data_readers/ CAsciiReader
	//{
	//
	// };
	bool write(string filename)
		{
		std::ofstream stream(filename.c_str());
		if(stream.is_open())
			{
			stream<<"# "<<m_comments<<endl;
			TDataMap::iterator it;
			for ( it=m_data.begin() ; it != m_data.end(); it++ )
				{
				stream << (*it).first << "\t"; 
				for(size_t i=0;i<(*it).second.size();i++)
					stream<<(*it).second[i]<<"\t";
				stream<<endl;
				}

			}else{
				stream.close();
				return false;
			}

		stream.close();
		return true;
		}
	void save(){
		
		if(!write(m_filename))
			cout<<"WARNING: cannot open file :"<<m_filename<<endl;
		else
			if(!write("dump.log"))
				{
				cout<<"ERROR: cannot open log file for writing:"<<"dump.log"<<endl;
				}
			
		};
inline	void SetComment(std::string comment){m_comments=comment;};
private:
	bool m_write_on_save;
	bool m_append;
	std::string m_comments;
	};

typedef std::vector< std::vector<float> > TLogData;

class CLoggerAm{

	typedef std::map<int, TLogData> TLogDataMap;
public:
	CLoggerAm(std::string file, bool update=false):m_filename(file), m_update(update)
		{
		  load();
		}
	~CLoggerAm(){save();}
	void flush(){save();};
	void SetComment(std::string msg){m_comment=msg;};
	bool insert(const int isnap, TLogData &data)
		{
		bool ret=!(m_data.find(isnap)!=m_data.end());
		for(size_t i=0;i<data.size();i++)
			m_data[isnap].push_back(data[i]);
		flush();
		return ret;
		}
	void save()
		{
		if(!write(m_filename))
			cout<<"WARNING: cannot open file :"<<m_filename<<endl;
		else
			if(!write("dump.log"))
				{
				cout<<"ERROR: cannot open log file for writing:"<<"dump.log"<<endl;
				}
		}

	bool load()
		{
		std::ifstream stream(m_filename.c_str());
		if(stream.is_open())
			{
			
			int isnap;
			char ch;
			std::string oneline;
			std::getline(stream,oneline, '\n');//skip comments
			//while(std::getline(stream,oneline, '\n'))
				{
				
				for(;;)
					{
					if(!std::getline(stream,oneline, '\n'))break;
					if(count_nonblanks(oneline)==0)continue;
					istringstream in(oneline);
					in>>ch>>isnap;
					
					
					//fill one profile
					for(;;)
						{
						std::vector<float> vec;
						std::getline(stream,oneline, '\n');
						if(count_nonblanks(oneline))
							{
							istringstream in(oneline);
							std::copy(std::istream_iterator<float>(in),std::istream_iterator<float>(), std::back_inserter(vec)); 
							m_data[isnap].push_back(vec);
							}else break;
						}
					}
					
				}

			}else{
				stream.close();
				return false;
			}

		stream.close();
		return true;
		}
private:
	bool write(string filename)
		{
		std::ofstream stream(filename.c_str());
		if(stream.is_open())
			{
			stream<<"# "<<m_comment<<endl;
			TLogDataMap::iterator it;
			for ( it=m_data.begin() ; it != m_data.end(); it++ )
				{
				stream <<"# "<< (*it).first << "\n"; 
				TLogData::iterator itdata=(*it).second.begin();
				size_t nj=(*it).second.size();
				size_t ni=(*it).second.begin()->size();
				for(size_t i=0;i<ni;i++)
					{
					for(size_t j=0;j<nj;j++)
						{
						stream<<(*it).second[j][i]<<"\t";
						}
					stream<<endl;
					}
				stream<<"\n\n"<<endl;
				}

			}else{
				stream.close();
				return false;
			}
		
		stream.close();
		return true;
		}
	//CLogger m_log;
	string m_filename;
	string m_comment;
	bool m_update;
	TLogDataMap m_data;
	};


#endif


