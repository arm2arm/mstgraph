#ifndef _APDATA_
#define _APDATA_

struct TApData{
	float AP[2];
	const float getAp(){return AP[1]/AP[0];};
	TApData(float A,float P){AP[0]=A;AP[1]=P;};
	TApData(){AP[0]=0;AP[1]=0;};
	};

#endif