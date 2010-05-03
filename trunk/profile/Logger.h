#ifndef _LOGGER_
#define _LOGGER_
#include "loader.h"
#include "utils.h"
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
	CLogger(std::string fout="dump.log", bool append=false,int nfields=1):m_append(append),CAsciiReader(fout, nfields,append ){


				    
		};
	~CLogger(){save();};
	////////// Populate ////////

	bool insert(const int isnap, const dynvector &d)
		{
		bool retval=is_done(isnap);
		m_data[isnap]=d;
		return retval;
		};
        bool is_done(unsigned int isnap)
	{return m_append&&(m_data.find(isnap)!=m_data.end());}
	//////////// IO ////////
	/*	void load()//moved to data_readers/ CAsciiReader
		{
		 ReadAsciiFile(); 
		 };*/
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
	void SetComment(std::string comment){m_comments=comment;};
private:
	bool m_append;
	std::string m_comments;
	};


#endif


