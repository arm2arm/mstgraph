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


enum eOVERPLOT{NO_ERASE, ERASE};
enum eCOLORS{BLACK_WHITE, BLUE_WHITE, GRN_RED_BLU_WHT, RED_TEMP, BGRY,STD_GAMMA, PRISM, RED_PURPLE, GREEN_WHITE_LIN, GREEN_WHITE_EXP,RAINBOW=13,HAZE=16, MAC_STYLE=25, RAINBOW_2=49  };

class CHOP{
public:
	CHOP(MyFloat alpha, int mingrp):m_alpha(alpha), m_mingrp(mingrp)
		{
		disp.assign(640*2,640,"Galaxy points",0);
		m_ctab=GREEN_WHITE_EXP;
		doRun();
		};
	unsigned int get_seeds(){return m_Nseeds;};
	/*
	Finally try to write the results.
	*/
	void write_catalogues(std::string hopfile)
		{

		}
private:
	CImgDisplay disp;//(640*2,640,"Galaxy points",0);
	CImg<unsigned char> A;//(disp.width()/2,disp.height(),1,3,255);
	CImg<unsigned char> B;//(disp.width()/2,disp.height(),1,3,255);
	int m_ctab;// fixed colortable for overploting;default is green table;
	void plotbyIndex(std::vector<int> ind)
		{
		std::vector<float> xyc;
		xyc.resize(ind.size()*5);
		float xrange[]={-5,5};
		int i=0;	
		for(unsigned int id=0;id<ind.size();id++)
			{
			i=ind[id];
			xyc[id]=Part[i].Pos[0];xyc[id+ind.size()]=Part[i].Pos[1];
			xyc[id+ind.size()*2]=Part[i].Est;
			xyc[id+ind.size()*3]=Part[i].Rho;
			xyc[id+ind.size()*4]=Part[i].Hsml/(xrange[1]-xrange[0])*640;
			}
		
		DrawVec(xyc, 5, xrange, xrange);
		}
	void oplotbyIndex(std::vector<int> ind, int cpal=GREEN_WHITE_EXP)
		{
		std::vector<float> xyc;
		xyc.resize(ind.size()*5);
		float xrange[]={-5,5};
		int i=0;	
		m_ctab=cpal;// set the global color table 
		for(unsigned int id=0;id<ind.size();id++)
			{
			i=ind[id];
			xyc[id]=Part[i].Pos[0];xyc[id+ind.size()]=Part[i].Pos[1];
			xyc[id+ind.size()*2]=256;
			xyc[id+ind.size()*3]=256;
			xyc[id+ind.size()*4]=Part[i].Hsml/(xrange[1]-xrange[0])*640;
			}

		DrawVec(xyc, 5, xrange, xrange, NO_ERASE);
		}
	int DrawVec(std::vector<float> xyc,const uint ncomp,  const float xrange[2], const float yrange[2],eOVERPLOT erase=ERASE) 
		{
		const unsigned int s_nb = xyc.size()/ncomp;
		float s_xmin = float(xrange[0]), s_xmax = float(xrange[1]);
				
		CImg<> samples(s_nb,2) ;
		CImg<float>color(s_nb),color1(s_nb);
		float *hsml=new float[s_nb];
		cimg_forX(samples,i) {			
			samples(i,0) = (float)(xyc[i]);
			samples(i,1) = (float)(xyc[i+s_nb]);
			color(i,0)   = (float)(xyc[i+s_nb*2]);
			color1(i,0)  = (float)(xyc[i+s_nb*3]);
			hsml[i]  = (float)(xyc[i+s_nb*4]);
			}

		disp.assign(640*2,640,"Galaxy points",0);
		CImg<unsigned char> Rhocol256, Estcol256;
		if(ERASE){
			color=color.log10();
			float meanCol=color.mean();
			float sc=255/(color.max()-meanCol);		
			color=(color-meanCol)*sc;

			color1=color1.log10();
			meanCol=color1.mean();
			color1=(color1-meanCol)/(color1.max()-meanCol)*255;
			cimg_forX(color,i){
				if(color(i)<2)color(i)=2;
				if(color1(i)<2)color1(i)=2;	
				}

			Rhocol256=color.normalize(0,255);
			Estcol256=color1.normalize(0,255);
		
			}else{
				Rhocol256=color;
				Estcol256=color1;
			}

		int ctab=13, i=0;
		unsigned int xc=(unsigned int)(disp.width()/2-5), yc=(unsigned int)(disp.height()*0.1-5);
		int bx1=5, by1=7;int bx2=xc, by2=yc;
		float btrans=0.5f;

		const unsigned char black[] = { 0,0,0 }, gray[] = { 228,228,228 };	
		const int maxLUT=CImg<int>::get_MAX_NUM_LUT();
		for (;!disp.is_closed() && !disp.is_keyQ() && !disp.is_keyESC();) 
			{    
			if(erase==ERASE)
				{
				A.assign(disp.width()/2,disp.height(),1,3,255);
				B.assign(disp.width()/2,disp.height(),1,3,255);
				B.fill(0);
				const unsigned int button = disp.button();
				if(button&1)
					ctab=(ctab+1)%maxLUT;
				if(button&2)
					ctab=(ctab-1)%maxLUT;
				}else
				{
				ctab=m_ctab;	
				};
			// Display points.							
			A.draw_galaxypoints(samples,Rhocol256,xrange, yrange, ctab);			
				//display(disp.resize(false).wait(1));
			CImg<unsigned char> bar=CImg<unsigned char>::IDL_LUT256(ctab);
			A.draw_rectangle(bx1-1,by1-1,bx2+1,by2+1,gray,btrans).draw_rectangle(bx1-1,by1-1,bx2+1,by2+1,black,1,~0U);
			for(uint i=0;i<256-1;i++)
				{
				const unsigned char col[]={bar(0,i,0),bar(0,i,1),bar(0,i,2)};
				A.draw_rectangle((int)(bx1+i*bx2/256.),by1,(int)(bx1+(i+1)*bx2/256.-1),by2,col,btrans);			
				}
			
			 //B.draw_galaxysprites(samples,Estcol256,xrange, yrange, ctab, hsml);
			B.draw_galaxypoints(samples,Estcol256,xrange, yrange, ctab, hsml);
			(A,B).display(disp);
			

			}
		A.save_jpeg("c:/arm2arm/DATA/testA.jpg");
		B.save_jpeg("c:/arm2arm/DATA/testB.jpg");
		return 0;
		}

	/////////////////////////////////////////
	typedef struct tagASTR{
		tagASTR():mmGrp(0),grpID(-1){};
		std::vector<uint> setB;
		int grpID;
		uint mmGrp;//most massive partner group
		} TMyNgb;

	vector<TMyNgb> NPart;// Particles linked to Groups 
	
	void setup_sets(std::vector<int> &isort=isortEst, eDOSORT mode=BY_EST)
		{
		uint loc_grpID=1, i=0;
		vector<int> iplot, iplot1;
		NPart.resize(All.NumPart);
		cout<<"Max:"<<Part[isort[0]].get(mode)<<endl;
		cout<<"Min:"<<Part[isort[All.NumPart-1]].get(mode)<<endl;
		for(uint ip=0;ip<NPart.size();ip++)
			{
			i=isort[ip];
			for(int ingb=0;ingb<3;ingb++)
				{
				if( Part[i].get(mode) < Part[ Part[i].getNGB(ingb) ].get(mode) )
					NPart[i].setB.push_back(Part[i].getNGB(ingb));
				}
				iplot.push_back(i);
			if(NPart[i].setB.size()==0)
				{
				NPart[i].grpID=loc_grpID++;				
				iplot1.push_back(i);
				}
			}
		m_Nseeds=count_seeds();
		std::reverse(iplot.begin(), iplot.end());
		std::reverse(iplot1.begin(), iplot1.end());
		plotbyIndex(iplot);
		iplot.clear();
		plotbyIndex(iplot1);
		iplot1.clear();
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

	void assign_12particles(void )
		{
		uint idx, iEst1, iEst2;
		id_4_passB.clear();
		for(uint i=0;i<NPart.size();i++)
			{
			idx=isortEst[i];
			if(NPart[idx].setB.size()==1)// here we need to assign particles to the mother group.
				{
				iEst1=NPart[idx].setB[0];
				NPart[idx].grpID=NPart[iEst1].grpID;
				//cout<<NPart[iEst1].grpID<<" "<<endl;
				}else
				{		
				if(NPart[idx].setB.size()==2)
					id_4_passB.push_back(idx);
				 }
			}
		oplotbyIndex(id_4_passB);
		for(uint i=0;i<id_4_passB.size();i++)
			{
			idx=id_4_passB[i];

			iEst1=NPart[idx].setB[0];
			iEst2=NPart[idx].setB[1];
			/*cout<<Part[iEst1].Est<<" "<<Part[iEst2].Est<<" "<<endl;
			cout<<NPart[iEst1].grpID<<" "<<NPart[iEst2].grpID<<" "<<endl;*/
			if(NPart[iEst1].grpID == NPart[iEst2].grpID)
				NPart[idx].grpID=NPart[iEst1].grpID;
			else
				id_4_passC.push_back(idx);
			}
		oplotbyIndex(id_4_passC, RED_TEMP);

		}
	void assign_2particles(void )
		{
		//uint idx, iEst1, iEst2;
		}
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


	void populate_structures()
		{

		uint i, nBad=0;
		cats.clear();
		cats.resize(m_Nseeds);
		for(i=0;i<NPart.size();i++)
			{
			if(NPart[i].grpID==0)
				{
				nBad++;
				}else{
					cats[NPart[i].grpID-1].id.push_back(i);
					cats[NPart[i].grpID-1].m_Np=cats[NPart[i].grpID-1].id.size();
				}


			}
		std::sort(cats.begin(),cats.end());
		std::reverse(cats.begin(),cats.end());
		SaveCats("c:/arm2arm/DATA/test.idx");
		cout<<"bad points: "<<nBad<<endl;
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
	/*********LOCAL variables */
	vector<TStr> cats;//Catalogue
	int m_mingrp;//minimum number of particles in the group 
	MyFloat m_alpha;//connectivity parameter.
	unsigned int m_Nseeds;//number of initial seeds for the groups
	vector<int> id_4_passB,id_4_passC;
	//////////// MAIN ALGORITHM DONE HERE //////////////////
	void doRun()
		{

		// pass A
			{
			scoped_timer timethis("#CHOP:> pass A: setup_sets\t"); 
			setup_sets(isortRho,BY_RHO);// lets fill the sets: setA and setB per particle
			}
			// pass B
			{
			scoped_timer timethis("#CHOP:> pass B: assign_12particles\t"); 
			assign_12particles();// lets fill the sets: setA and setB per particle
			}
			// pass C
			{
			scoped_timer timethis("#CHOP:> pass C: assign_2particles\t"); 
			//assign_2particles();// lets fill the sets: setA and setB per particle
			}
			//populate_structures();
		}


	};//end of CHOP class


#endif