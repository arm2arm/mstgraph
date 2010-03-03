#include "HOP.h"

/*
void CHOP::plotbyIndex(std::vector<int> ind){
	std::vector<float> xyc;
	xyc.resize(ind.size()*5);
	float xrange[]={-5,5};
	int i=0;	
	for(unsigned int id=0;id<ind.size();id++){
		i=ind[id];
		xyc[id]=Part[i].Pos[0];xyc[id+ind.size()]=Part[i].Pos[1];
		xyc[id+ind.size()*2]=Part[i].Est;
		xyc[id+ind.size()*3]=Part[i].Rho;
		xyc[id+ind.size()*4]=Part[i].Hsml/(xrange[1]-xrange[0])*640;
		}

	DrawVec(xyc, 5, xrange, xrange);
	}
void CHOP::oplotbyIndex(std::vector<int> ind, int cpal){
	std::vector<float> xyc;
	xyc.resize(ind.size()*5);
	float xrange[]={-5,5};
	int i=0;	
	m_ctab=cpal;// set the global color table 
	for(unsigned int id=0;id<ind.size();id++){
		i=ind[id];
		xyc[id]=Part[i].Pos[0];xyc[id+ind.size()]=Part[i].Pos[1];
		xyc[id+ind.size()*2]=256;
		xyc[id+ind.size()*3]=256;
		xyc[id+ind.size()*4]=Part[i].Hsml/(xrange[1]-xrange[0])*640;
		}

	DrawVec(xyc, 5, xrange, xrange, NO_ERASE);
	}
int CHOP::DrawVec(std::vector<float> xyc,const uint ncomp,  const float xrange[2], const float yrange[2],eOVERPLOT erase) 
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
			for(uint i=0;i<256-1;i++){
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


*/

