#ifndef _MYREADERS_
#define _MYREADERS_

#ifndef M_PI
#define M_PI 3.141592653589793238462643383
#endif
#define ND 6

typedef int int4byte;
#include <vector>

#define real double
#define real1 float

void read_ic12(const char *fname);
void read_ic12_est(const char *fname);
void read_ic12_rho(const char *fname);
void write_block(std::vector<float> vec,const char *fname, char *TAG_NAME);
void fill_isort_vecs();

void SwapEndian(void* addr, int* pattern);
void   endrun(int);
void allocate_memory(void);
double drand48();
/* this struct contains mostly code parameters read from the parameter file */
extern struct global_data_all_processes
{
    /* Code options */
    int    TypeOfSmoothing;
    int    KernelBiasCorrection;
    int    AnisotropicKernel;
    int    Dimensions;
    int    VolCorr;
    int    TypeOfKernel;
    double SpatialScale;
    int    PartBoundary;
    int    NodeSplittingCriterion;
    int    CubicCells;
    double Anisotropy;
    int    DesNumNgb;
    int    NumBucket;
    int    DesNumNgbA;
    int    NumBucketA;
    /* Cosmology */
//    double BoxSize, BoxHalf;
    int    MedianSplittingOn;
    float   hs[ND],hsv;
    double  MassTable[6];
    /* File options */
    int   ICFormat;
    double OmegaBaryon;
    double Omega0;
    char    InitCondFile[100],SnapshotFileBase[100];
    char   InputDir[100],
	InitFileBase[100];
    /* Some other global parameters */
    int   MaxPart,TypeListOn;
    int flag_swap;
    int order_flag;
    int PeriodicBoundaryOn;
	int ONLY_TYPE;
	int NumPart;
#ifdef PERIODIC
    double boxh[ND];
#endif

} All;


extern std::vector<int> npart,npartc; 

static enum eDOSORT{BY_RHO, BY_EST, BY_POS};
extern struct particle_data
{	
    int  ID,id,Type;//,NumNgb;           /* unique particle identifier */
    float     Pos[ND];       // particle position
#ifdef DIM3
    float     Vel[3];       // particle velocity
#endif
    float     Mass;         /* particle mass */
    float     Est;
	float     Rho;
	float     Hsml;
	int  *pNGB;	
	int *pNGBR;
//	float get(eDOSORT mode=BY_EST);
//	int getNGB(int i,eDOSORT mode=BY_EST);
inline float get(eDOSORT mode=BY_EST){ if(mode==BY_EST)return Est;return Rho;};
inline int   getNGB(int i,eDOSORT mode=BY_EST){if(mode==BY_EST)return pNGB[i];return pNGBR[i];};

} *P,*P_data,*Part;//,*P1;


extern std::vector<int> isortRho,isortEst, isortPos;

/* Header for the standard file format. */
extern struct io_header_1
{
    int4byte npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int4byte flag_sfr;
    int4byte flag_feedback;
    int4byte npartTotal[6];
    int4byte flag_cooling;
    int4byte num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    /* extra flags other than gadget*/
    int4byte flag_id;          /* if IDs needed in output */
    int4byte flag_dim;         /* no of dimensions */
    /* signifies that output file contains estimated density */
    int4byte flag_density;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8-3*4];  /* fills to 256 Bytes */
} header1;

#endif