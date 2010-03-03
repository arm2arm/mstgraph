
#include "Render.h"
#include <cmath>
#include <algorithm>
#ifdef WITH_CIMG
CRender::CRender(void)
	{
	}

CRender::~CRender(void)
	{
	}
void DoJob()
	{


	}

void CRender::setup_cic_image()
	{
	img.assign(512,512,1,3,0);
	if(m_visu)
		{
		
		CImg<float> ColorTable=CImg<float>::rainbow_LUT256();
		int i;
		int w=img.width();
		int h=img.height();
		MyFloat minX=1e10, maxX=-1e10, Cent[2]={0,0};
		int shape=1,profile=0;
		int x, y;
		MyFloat xf, yf, Rad=6, Rad2=Rad*Rad, r;
		int hsml;//Hsml  scaled to Image pixels
		// Define the color of the gradient
		CImg<unsigned char> col(3);
		const unsigned char col1[3] = { 0,0,255 }, white[3] = { 255,255,255 };
		CImgDisplay disp(w*2, w,"Click and drag to create color gradient",0);
		disp.set_title("Random color gradient").show().flush();
		for(i=1;i<=All.NumPart;i++)
			{
			Cent[0]+=P[i].Pos[0];
			Cent[1]+=P[i].Pos[1];
			minX=std::min<MyFloat>(minX,P[i].Pos[0]);
			maxX=std::max<MyFloat>(maxX,P[i].Pos[0]);
			}
		Cent[0]/=MyFloat(All.NumPart);
		Cent[1]/=MyFloat(All.NumPart);
		minX=-Rad;maxX=Rad;
		CImg<unsigned char> visu(img);
		int ixc, iyc, ip;
		MyFloat sphW=0;
		MyFloat miRho=log10(P[isortRho[All.NumPart-1]].Rho),maRho=log10(P[isortRho[0]].Rho);
		miRho=-7;
		MyFloat scaleRho=256.0/(maRho-miRho), Rho;
		MyFloat alpha=0.5, meanRho=(maRho-miRho)/2;
		for(i=All.NumPart;i>=1;i--)
			{

			if( drand48()<0.01)
				{
				ip=isortRho[i-1];
				xf=(P[ip].Pos[0]-Cent[0]);
				yf=(P[ip].Pos[1]-Cent[1]);
				x=(int)((xf+Rad)/(2*Rad)*w);
				y=(int)((yf+Rad)/(2*Rad)*h);
				if(x>0 && x<w && y<w && y>0)
					{
					hsml=int(P[ip].Hsml/(Rad)*h);
					Rho=(log10(P[ip].Rho)-miRho)*scaleRho;

					//hsml=10;
					//sphW=0;
					//img.draw_gradient(x,y,x+hsml,y+hsml,col1,0,shape,profile,.4f);
					for(int ix=-hsml;ix<hsml+1;ix++)
						for(int iy=-hsml;iy<hsml+1;iy++)
							{
							ixc=x+ix;
							iyc=y+iy;
							r=std::sqrt(MyFloat(ix*ix+iy*iy));
							sphW=Rho*Wsph<MyFloat>(r, hsml);

							if(ixc<0)ixc=0;
							if(ixc>w)ixc=w-1;
							if(iyc<0)iyc=0;
							if(iyc>w)iyc=w-1;
							//GetColor(sphW);
							//if(sphW > meanRho)
								{							
								visu(ixc,iyc)=(int)((sphW+visu(ixc,iyc))*0.5);//Red
								visu(ixc,iyc,0,1)=visu(ixc,iyc);
								visu(ixc,iyc,2)=visu(ixc,iyc);

								}
								
							}
						
						if(i % 100 ==0)cout<<"i="<<i<<" "<<Rho<<"  Img("<<x<<","<<y<<")="<<img(x,y)<<"\t\t\r";
						/*alpha=Rho;
						cimg_forXYZC(img,cx,cy,cz,ck) {
							unsigned char pix=(unsigned char)((1 - alpha)*visu(cx,cy,cz,ck) + alpha*img(cx,cy,cz,ck));
						 img(cx,cy,cz,ck) = pix ;
						}
						//img_forXYZC(visu,cx,cy,cz,ck) {visu(cx,cy,cz,ck)=0;}
						//img.draw_gaussian(x,y,0.5,col1, 0.5);
						//img.draw_gradient(x,y,x+hsml,y+hsml,col1,0,shape,profile,.04f);
						*/
						(visu,img).display(disp);
					}
					//cout<<xf<<" "<<yf<<" "<<x<<" "<<y<<"\n";
					/*visu = img.normalize(0,255);
					visu.draw_text(10,10,"%.1ffps",col1,0,1,11,disp.frames_per_second()).display(disp);
					if (disp.is_resized()) disp.resize();*/
					

				}


			}
		cout<<endl;
		
		
		//img=img.normalize(img.mean(),maRho);
		
//		img/=img.mean();

		const CImg<float> imnorm=img.normalize(0,255);

		imnorm.display();
		

		}
	}

void CRender::PlotData(std::vector<MyFloat>  &vec)
	{
	// Create plot data
	int Nbins=256, inx;
	unsigned int plot_type=1;
	unsigned int vertex_type=1;

	CImg<double> values(Nbins,Nbins);
	MyFloat mi=1.0e-8;//*std::min_element(vec.begin(),vec.end());
	MyFloat ma=*std::max_element(vec.begin(),vec.end());
	for(int i=0;i<Nbins;i++)values(i)=i;  
	for(unsigned int iv=0;iv<vec.size();iv++)
		{
		inx=(int)((vec[iv]-mi)/(ma-mi)*Nbins);
		if(inx >-1 && inx<(Nbins-1))
			{
			values(0,inx) +=1;

			}
		}
	values.spectrum();
	values.display_graph(0,3);
	}
#endif
///////////////////////////////
