#ifndef _MY_HOP_
#define _MY_HOP_
/*
This is the main algorithm to find the FOF like groups in the 6D phase space.
this requre readers.h for the particle data.
*/
#include <cstring>
#include <vector>
#include <map>
#include "functions.h"
#include "readers.h"
#include "global_typedefs.h"
#include <vector>
#include "cats.h"
using std::vector;
using std::cout;
using std::endl;


enum eOVERPLOT{NO_ERASE, ERASE};
enum eCOLORS{BLACK_WHITE, BLUE_WHITE, GRN_RED_BLU_WHT, RED_TEMP, BGRY,STD_GAMMA, PRISM, RED_PURPLE, GREEN_WHITE_LIN, GREEN_WHITE_EXP,RAINBOW=13,HAZE=16, MAC_STYLE=25, RAINBOW_2=49  };

class CHOP{
public:
	CHOP(MyFloat alpha, int mingrp):m_alpha(alpha), m_mingrp(mingrp)
		{
		//disp.assign(640*2,640,"Galaxy points",0);
		//m_ctab=GREEN_WHITE_EXP;
		doRun();
		};
	unsigned int get_seeds(){return m_Nseeds;};
	/*
	Finally try to write the results.
	*/
	void write_catalogues(std::string hopfile)
		{

		};
private:
	//CImgDisplay disp;//(640*2,640,"Galaxy points",0);
	//CImg<unsigned char> A;//(disp.width()/2,disp.height(),1,3,255);
	//CImg<unsigned char> B;//(disp.width()/2,disp.height(),1,3,255);
	//int m_ctab;// fixed colortable for overploting;default is green table;
	//void plotbyIndex(std::vector<int> ind);		
	//void oplotbyIndex(std::vector<int> ind, int cpal=GREEN_WHITE_EXP);		
	//int DrawVec(std::vector<float> xyc,const uint ncomp,  const float xrange[2], const float yrange[2],eOVERPLOT erase=ERASE); 

	/////////////////////////////////////////
	typedef struct tagASTR{
		tagASTR():mmGrp(0),grpID(-1),i(0){};
		std::vector<uint> setB;
		int grpID;
		int i;
		uint mmGrp;//most massive partner group
		} TMyNgb;

	vector<TMyNgb> NPart;// Particles linked to Groups 
	vector<int> iplot, iplot1, set0,set1, set2, set3;
	void setup_sets(std::vector<int> &isort=isortEst, eDOSORT mode=BY_EST)
		{
		uint loc_grpID=1, i=0, j=0;
		iplot.clear();		
		iplot1.clear();
		NPart.resize(All.NumPart);
		cout<<"Max:"<<Part[isort[0]].get(mode)<<endl;
		cout<<"Min:"<<Part[isort[All.NumPart-1]].get(mode)<<endl;
		for(uint ip=0;ip<NPart.size();ip++)
			{
			i=isort[ip];
			for(int ingb=0;ingb<2;ingb++)
				{
				j=Part[i].getNGB(ingb);
				if( Part[i].get(mode) < Part[j].get(mode) )
					NPart[i].setB.push_back(j);
					NPart[i].i=ip;
				}
				iplot.push_back(i);
			if(NPart[i].setB.size()==0)
				{
				NPart[i].grpID=loc_grpID++;				
				iplot1.push_back(i);
				set0.push_back(i);
				}
			if(NPart[i].setB.size()==1)
				set1.push_back(i);
			if(NPart[i].setB.size()==2)
				set2.push_back(i);
			if(NPart[i].setB.size()==3)
				set3.push_back(i);
			}
		m_Nseeds=count_seeds();
		std::reverse(iplot.begin(), iplot.end());
		std::reverse(iplot1.begin(), iplot1.end());
		//plotbyIndex(iplot);
		//iplot.clear();
		//plotbyIndex(iplot1);
		//iplot1.clear();
		printStat();
		}
void printStat()
	{
	cout<<"#############"<<endl;
	cout<<"setB=0 :"<<set0.size()<<endl;
	cout<<"setB=1 :"<<set1.size()<<endl;
	cout<<"setB=2 :"<<set2.size()<<endl;
	cout<<"setB=3 :"<<set3.size()<<endl;
	cout<<"#############"<<endl;
	}
	unsigned int count_seeds()
		{
		unsigned int i, counter=0;
		vector<uint> seed_vector;
		for(i=0;i<NPart.size();i++)
			{
			if(NPart[i].setB.size()==0)
				{
				 counter++;
				}
			}		
		return counter;
		}

	void assign_12particles(std::vector<int> &isort=isortEst, eDOSORT mode=BY_EST)
		{
		uint i1, i2,i;
		id_4_passB.clear();
		for(uint ip=0;ip<NPart.size();ip++)
			{
			i=isort[ip];
			if(NPart[i].setB.size()==1)// here we need to assign particles to the mother group.
				{
				i1=NPart[i].setB[0];
				if(NPart[i1].grpID==-1)
					assert(0);
				NPart[i].grpID=NPart[i1].grpID;				
				}else
				{		
				if(NPart[i].setB.size()==2)
					id_4_passB.push_back(i);// we have two particles in the set B eg: to which particle to attach?
				else
					if(NPart[i].setB.size()==3)
						id_4_passC.push_back(i);// an alternative attachement 
					else
						assert(13);
				 }
			}
		
		for(uint ip=0;ip<id_4_passB.size();ip++)
			{
			i=id_4_passB[ip];

			i1=NPart[i].setB[0];
			i2=NPart[i].setB[1];
			//if(NPart[i1].grpID ==-1 && NPart[i2].grpID ==-1)assert(0);// something is goes wrong
			if(NPart[i1].grpID == NPart[i2].grpID )
				NPart[i].grpID=NPart[i1].grpID;
			else
				{
					//if(Part[i1].get() > Part[i2].get())

						
				}
			}
		
		}
	void assign_2particles(void )
		{
		uint i,i1,i2,i3;
		for(uint ip=0;ip<id_4_passC.size();ip++)
			{
			i=id_4_passC[ip];

			i1=NPart[i].setB[0];
			i2=NPart[i].setB[1];
			i3=NPart[i].setB[2];
//		cout<<i<<endl;
			if(NPart[i1].grpID == NPart[i2].grpID && NPart[i3].grpID == NPart[i2].grpID && NPart[i1].grpID!=-1)
				NPart[i].grpID=NPart[i1].grpID;
			else
				id_4_passD.push_back(i);
			}
		}

	CCatalog c;
	void populate_structures()
		{

		uint i, nBad=0;
		c.cats.clear();
		c.cats.resize(m_Nseeds);
		for(i=0;i<NPart.size();i++)
			{
			if(NPart[i].grpID==-1)
				{
				nBad++;
				}else{
					c.cats[NPart[i].grpID-1].id.push_back(i);
					c.cats[NPart[i].grpID-1].m_Np=c.cats[NPart[i].grpID-1].id.size();
				}


			}
		c.sort();
		cout<<"bad points: "<<nBad<<endl;
		c.SaveCats("c:/arm2arm/DATA/test.idx");

		}
	/*********LOCAL variables */
	int m_mingrp;//minimum number of particles in the group 
	MyFloat m_alpha;//connectivity parameter.
	unsigned int m_Nseeds;//number of initial seeds for the groups
	vector<int> id_4_passB,id_4_passC,id_4_passD;
	//////////// MAIN ALGORITHM DONE HERE //////////////////
	void doRun()
		{

		// pass A
			{
			scoped_timer timethis("#CHOP:> pass A: setup_sets\t"); 
			setup_sets(isortEst,BY_EST);// lets fill the sets: setA and setB per particle
			}
			// pass B
			{
			scoped_timer timethis("#CHOP:> pass B: assign_12particles\t"); 
			assign_12particles(isortEst,BY_EST);// lets fill the sets: setA and setB per particle
			}
			// pass C
			{
			scoped_timer timethis("#CHOP:> pass C: assign_2particles\t"); 
	//		assign_2particles();// lets fill the sets: setA and setB per particle
			}
			populate_structures();
		}


	};//end of CHOP class


#endif