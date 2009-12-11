//=======================================================================
// Copyright 2009 Arman Khalatyan, OAMP/LAM  
// distributed under GPL license
//=======================================================================
#include "mst_graph.h"
#include "kdtree2.hpp"
#include "functions.h"
#include "Render.h"
#include "HOP.h"
#include "GetEst.h"
/////////////////////////////
//#define ND_GROUPS 2
#define LOAD_DEBUG 0
#define SMALL_TEST 0
#define MIN_NGRP 50
#define MAXNGB 64.0 
std::string base;//this is used for outputs
#define RAND_FRACTION 1.2
///////////////////////////////
//
// define, for convenience a 2d array of floats. 
//  
typedef boost::multi_array<MyFloat,2> array2dfloat;

///////////////////////////////
#include "readers.h"
void read_gadget(string fname)
	{
	All.ONLY_TYPE=4;
	read_ic12(fname.c_str(), false);
/*	MyFloat COM[]={0.0,0.0,0.0}, totEst=0.0;
	for(int i=0;i<All.NumPart;i++)
		{
		COM[0]+=Part[i].Pos[0]*Part[i].Est;		
		COM[1]+=Part[i].Pos[1]*Part[i].Est;		
		COM[2]+=Part[i].Pos[2]*Part[i].Est;		
		totEst+=Part[i].Est;
		}
	COM[0]/=totEst;
	COM[1]/=totEst;
	COM[2]/=totEst;
	for(int i=0;i<All.NumPart;i++)
		{
		Part[i].Pos[0]-=(float)COM[0];
		Part[i].Pos[1]-=(float)COM[1];
		Part[i].Pos[2]-=(float)COM[2];
		}
		*/
	}
///////////////////////////////////////////////////////////////////
//"snap_gal_sfr_0450.ascii"	
int main(int argc,char **argv) {

	// Read and check command line parameters.
	cimg_usage("Compute a HOP over the particles with given Est and Rho files");
	const char *file_ascii = cimg_option("-ascii_file",(char*)0,"Input in Ascii format");
	const char *file_gadget2;
	if( SMALL_TEST)
		file_gadget2= cimg_option("-gadget2_file","C:\\arm2arm\\DATA\\MODEL7\\MODELS\\MODEL7\\RUNG2\\SNAPS\\test_0450","Input in Gadget-2 format");
	else
		file_gadget2= cimg_option("-gadget2_file","C:\\arm2arm\\DATA\\MODEL7\\MODELS\\MODEL7\\RUNG2\\SNAPS\\snap_gal_0.25_sfr_0303","Input in Gadget-2 format");
	const char *file_base  = cimg_option("-base4files","C:/arm2arm/DATA/MODEL7/MST_GRAPH/","Base for input and output");
	const char *HOP_file  = cimg_option("-hopfile","hop_0450","HOP algorithm output,eg: groups and their IDs.");
	const bool visu     = cimg_option("-visu",true,"Visualization mode");
	const int shape  = cimg_option("-s",1,"shape [0,6]");
	const int profile  = cimg_option("-p",0,"profile [0,7]");

	////////////////////////////////////////////////////////////////////

	base=string(file_base);
	string hopfile=base+string("/")+string(HOP_file);
		
	read_gadget(file_gadget2);

	////////////////////////////////////////////////
	/// Smooth Est with 64 Ngb
	GetEst est_me;
	est_me.Run_SPHEst();
	exit(0);
	//SmoothSph();
	/////////////////////// GET setAB This is the HOP stuff based on Enbid Density  one can test also for RHO by SPH//////////////////////////////////////
	//exit(0);
	MyFloat alpha=0.05;
	CHOP *hop=new CHOP(alpha, MIN_NGRP);

	int seeds=hop->get_seeds();
	cout<<"We got "<<seeds<<" seeds out of "<< All.NumPart<<endl;


	hop->write_catalogues(hopfile);

	///////////////////////////////////////////////////////////////
	// Here we need to cut and store the catalogues
	///////////////////////////////////////////////////////


	////////////////////////////////
	return EXIT_SUCCESS;
	}
