#include <stdio.h>
#include<iostream>
#include<vector>
#include <algorithm>
#include "readers.h"
using namespace std;
////////////////
/* this struct contains mostly code parameters read from the parameter file */
struct global_data_all_processes All;

struct io_header_1 header1;

//struct param_data Var;

struct particle_data *P, *P_data, *Part;



/* The following structure holds all the information that is
 * stored for each particle.
 */

vector<int> npart, npartc;
vector<int> isortRho, isortEst, isortPos;

////////////////

size_t my_fwrite1(void *ptr, int* pattern, size_t nmemb, int flag_swap, FILE *stream) {
    size_t nwritten, size, i;
    //    int i;
    void* addr;

    for (i = 0, size = 0; pattern[i] > 0; i += 2)
        size += pattern[i] * pattern[i + 1];

    if (flag_swap == 1) {
        for (i = 0; i < nmemb; i++) {
            addr = (char *) ptr + i*size;
            SwapEndian(addr, pattern);
        }
    }


    if ((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb) {
        printf("I/O error (fwrite) on has occured.\n");
        fflush(stdout);
        endrun(777);
    }

    if (flag_swap == 1) {
        for (i = 0; i < nmemb; i++) {
            addr = (char *) ptr + i*size;
            SwapEndian(addr, pattern);
        }
    }

    return nwritten;
}

size_t my_fread1(void *ptr, int *pattern, size_t nmemb, int flag_swap, FILE *stream) {
    size_t nread, i;
    void *addr;
    int size;

    for (i = 0, size = 0; pattern[i] > 0; i += 2)
        size += pattern[i] * pattern[i + 1];


    if ((nread = fread(ptr, size, nmemb, stream)) != nmemb) {
        printf("I/O error (fread) has occured.\n");
        fflush(stdout);
        endrun(778);
    }


    if (flag_swap == 1) {
        for (i = 0; i < nmemb; i++) {
            addr = (char *) ptr + i*size;
            SwapEndian(addr, pattern);
        }
    }

    return nread;
}

void SwapEndian(void* addr, int* pattern) {
    int i, j = 0, k;
    char c;
    while (pattern[j] > 0) {
        for (k = 0; k < pattern[j + 1]; k++) {
            for (i = 0; i < pattern[j] / 2; i++) {
                c = *((char*) addr + i);
                *((char*) addr + i) = *((char*) addr + (pattern[j] - i - 1));
                *((char*) addr + (pattern[j] - i - 1)) = c;
            }
            addr = ((char *) addr) + pattern[j];
        }
        j += 2;
    }
}

void WriteOneBlock(FILE *file, string blname, char* pData, unsigned int datasize) {
    unsigned int blsize, idata;
    /*write block name*/
    blsize = 8;
    fwrite((char*) & blsize, 4, 1, file);
    fwrite(blname.c_str(), 4, 1, file);
    idata = 8 + datasize;
    fwrite((char*) & idata, 4, 1, file);
    fwrite((char*) & blsize, 4, 1, file);
    /////////////////////////////////
    /*Writing data*/
    blsize = datasize;
    fwrite((char*) & blsize, 4, 1, file);
    fwrite((char*) pData, datasize, 1, file);
    fwrite((char*) & blsize, 4, 1, file);

}

void write_block(vector<float> vec, const char *fname, char *TAG_NAME) {
    FILE *file;
    int nv = vec.size();

    if ((file = fopen(fname, "wb"))) {
        int blsize = sizeof (header1);
        WriteOneBlock(file, "HEAD", (char*) & header1, blsize);
        blsize = sizeof (float) * nv;
        WriteOneBlock(file, string("AGE "), (char*) (&vec[0]), blsize);
        fclose(file);
    }
}

/* GADGET 2 Format  */
void read_ic12(const char *fname, bool ngb_flag) {
#define GETBLKNAME fread(&blkname, sizeof(blkname), 1, fd);cout<<blkname.name[0]<<blkname.name[1]<<blkname.name[2]<<blkname.name[3]<<endl;
#define SKIP my_fread1(&blklen,patterns,1,All.flag_swap,fd);
    FILE *fd;
    int i, k, massflag, count;
    float dummy[3];

    struct {
        int blks;
        char name[4];
        int next;
        int blk2;
    } blkname;
    int pc, type;
    int patternd[] = {4, 3, 0}, patterns[] = {4, 1, 0}, patternh[] = {4, 6, 8, 8, 4, 10, 8, 4, 4, 3, 4, 21, 0};
    int4byte intdummy, blklen;
    int NumPart, ONLY_TYPE = All.ONLY_TYPE; // this will select only specific type of the particles

    All.flag_swap = 0;




    if ((fd = fopen(fname, "rb"))) {
        fprintf(stdout, "Reading GADGET format file: %s \n", fname);
        fflush(stdout);

        GETBLKNAME;
        SKIP;
        if (blklen != 256) {
            All.flag_swap = 1;
            SwapEndian(&blklen, patterns);
            if (blklen != 256) {
                printf("incorrect header format (2)\n");
                endrun(889);
            }
        }

        my_fread1(&header1, patternh, 1, All.flag_swap, fd);


        SKIP;
        if (blklen != 256) {
            printf("incorrect header format (2)\n");
            endrun(889);
        }

        //	All.BoxSize=header1.BoxSize;
        for (i = 0, massflag = 0; i < 6; i++) {
            if (header1.mass[i] == 0 && header1.npart[i] > 0)
                massflag = 1;
        }


        for (i = 0, NumPart = 0; i < 6; i++)
            NumPart += header1.npart[i];

        All.MaxPart = NumPart; // sets the maximum number of particles
        cout << "got total number of particles:" << NumPart << endl;
        allocate_memory();

        GETBLKNAME;
        SKIP;
        float *vec = new float[NumPart * 3];
        int ic, nread;
        nread = fread((char*) & vec[0], sizeof (float), NumPart * 3, fd);
        if (nread != NumPart * 3)
            endrun(778);
        for (i = 1, ic = 0; i <= NumPart; i++) {
            // my_fread1(&dummy[0],patternd,1,All.flag_swap,fd);
            for (k = 0; k < 3; k++)
                P[i].Pos[k] = vec[ic + k];
            ic += 3;
        }
        SKIP;

        GETBLKNAME;
        SKIP;
        for (i = 1; i <= NumPart; i++) {
            my_fread1(&dummy[0], patternd, 1, All.flag_swap, fd);
#ifdef DIM3
            for (k = 0; k < 3; k++)
                P[i].Vel[k] = dummy[k];
#else
            if (NDIM > 3)
                for (k = 0; k < (NDIM - 3); k++)
                    P[i].Pos[k + 3] = dummy[k];
#endif
        }
        SKIP;

        GETBLKNAME;
        SKIP;
        for (i = 1; i <= NumPart; i++) {
            my_fread1(&intdummy, patterns, 1, All.flag_swap, fd);
            P[i].ID = intdummy;
        }
        SKIP;


        if (massflag) {
            GETBLKNAME;
            SKIP;
        }
        for (type = 0, count = 1; type < 6; type++) {
            if (header1.mass[type] == 0 && header1.npart[type] > 0) {
                k = count;
                for (i = 1; i <= header1.npart[type]; i++) {
                    my_fread1(&dummy[0], patterns, 1, All.flag_swap, fd);
                    P[count++].Mass = dummy[0];

                }

            } else {
                for (i = 1; i <= header1.npart[type]; i++) {
                    P[count++].Mass = (float) header1.mass[type];
                }
            }
        }
        if (massflag)
            SKIP;

        fclose(fd);
        fprintf(stdout, "....done with reading.\n");
        fflush(stdout);

        /*Compress vectors to select only one type*/
        if (ONLY_TYPE > (-1)) {
            int Nonly = header1.npart[ONLY_TYPE], pshift = 0, mshift = 0;
            for (type = 0; type < ONLY_TYPE; type++) {
                pshift += header1.npart[type];

            }
            for (type = 0; type < 6; type++) {
                if (ONLY_TYPE != type)
                    header1.npart[type] = 0;
            }
            for (i = 1; i <= Nonly; i++) {
                memcpy(&P[i], &P[i + pshift], sizeof (struct particle_data));
                P[i].id = i - 1;
            }
            NumPart = Nonly;
            header1.npart[ONLY_TYPE] = NumPart;
            //cout<<P[1].Mass<<endl;
            //cout<<P[NumPart].Mass<<endl;
            for (i = 1; i <= NumPart; i++)
                if (P[i].Mass != P[1].Mass)
                    cout << i << " " << P[i].Mass << endl;
            for (i = 1; i <= NumPart; i++)
                P[i].Mass = 1;
        }
        ////////////////////////////////////////////

        npart.clear();
        npartc.clear();
        for (i = 0, k = 0; i < 6; i++) {
            npart.push_back(header1.npart[i]);
            npartc.push_back(k);
            k += npart[i];
        }
        /* set the particle types */
        for (type = 0, pc = 1; type < 6; type++)
            for (i = 0; i < npart[type]; i++)
                P[pc++].Type = type;


    } else {
        fprintf(stdout, "File %s not found.\n", fname);
        endrun(7);
    }

    for (i = 0; i < 6; i++)
        if (header1.npart[i] > 0) cout << "Type = " << i << " Particles = " << header1.npart[i] << endl;


    fprintf(stdout, "Total particles = %d\n", NumPart);
    All.NumPart = NumPart;

    char frho[1024], fest[1024];
    sprintf(&frho[0], "%s%s", fname, "_rho_4");
    sprintf(&fest[0], "%s%s", fname, "_ph4.est");
    if (ngb_flag) {
        //read_ic12_est(fest);
        read_ic12_rho(frho);
        //fill_isort_vecs();
    }
#undef SKIP    
#undef GETBLKNAME  
}
//End for GADGET 2 format



// sorter helper functions

class IntSequence {
private:
    int value;
public:
    // constructor

    IntSequence(int initialValue)
    : value(initialValue) {
    }

    // ''function call''

    int operator() () {
        return value++;
    }
};

template <class T>
struct setGtRho {
    T *that;

    setGtRho(T * f) : that(f) {
    }

    bool operator()(const unsigned int p1, const unsigned int p2) {
        return (that[p1].Rho > that[p2].Rho);
    }
};

template <class T>
struct setGtEst {
    T *that;

    setGtEst(T * f) : that(f) {
    }

    bool operator()(const unsigned int p1, const unsigned int p2) {
        return (that[p1].Est > that[p2].Est);

    }
};

template <class T>
struct setGtPos {
    T *that;

    setGtPos(T * f) : that(f) {
    }

    bool operator()(const unsigned int p1, const unsigned int p2) {
        return (
                that[p1].Pos[0] * that[p1].Pos[0] +
                that[p1].Pos[1] * that[p1].Pos[1] +
                that[p1].Pos[2] * that[p1].Pos[2]
                <
                that[p2].Pos[0] * that[p2].Pos[0] +
                that[p2].Pos[1] * that[p2].Pos[1] +
                that[p2].Pos[2] * that[p2].Pos[2]);
    }
};

/*make sort indexes for fast access*/
void fill_isort_vecs() {
    IntSequence seq(0); // integral sequence starting with 1

    std::generate_n<back_insert_iterator<vector<int> >,
            int, IntSequence&>(back_inserter(isortRho), // start
            All.NumPart, // number of elements
            seq); // generates values	}
    isortEst = isortRho;
    isortPos = isortRho;
    cout << "Fill the sort vectors" << endl;
    ///////////////// SORT by Rho ////////////
    cout << "NOT: Sorting by EST" << endl;
    //sort(isortEst.begin(),isortEst.end(),setGtEst<struct particle_data>(Part));
    cout << "Max:" << Part[isortEst[0]].Est << endl;
    cout << "Min:" << Part[isortEst[All.NumPart - 1]].Est << endl;
    cout << "NOT: Sorting by Rho" << endl;
    //sort(isortRho.begin(),isortRho.end(),setGtRho<struct particle_data>(P));
    cout << P[isortRho[0]].Rho << endl;
    cout << P[isortRho[All.NumPart - 1]].Rho << endl;

    float mean[] = {0, 0, 0, 0, 0, 0};
    int ii;
    cout << "Sorting by Ngb lists by Est" << endl;
    for (int i = 0; i < All.NumPart; i++) {
        //sort(&Part[i].pNGB[0],&Part[i].pNGB[All.DesNumNgb-1],setGtEst<struct particle_data>(Part));
        for (ii = 0; ii < 6; ii++)mean[ii] += Part[i].Pos[ii];

    }
    for (ii = 0; ii < 6; ii++)mean[ii] /= (float) All.NumPart;

    // remove mean

    for (int i = 0; i < All.NumPart; i++) {
        for (ii = 0; ii < 6; ii++)
            Part[i].Pos[ii] -= mean[ii];
        //		cout<<Part[i].Pos[0]<<" "<<Part[i].Pos[1]<<" "<<Part[i].Pos[2]<<endl;
    }

}

/* end the run */
void endrun(int ierr) {
    if (ierr) {
        fprintf(stdout, "endrun called with an error level of %d\n\n\n", ierr);
        exit(1);
    }
    exit(0);
}

/* This routine allocates memory for 
 * particle storage.
 */
void allocate_memory(void) {


    int bytes = 0, bytes_tot = 0;

    if (All.MaxPart > 0) {
        if (!(P_data = new particle_data[All.MaxPart])) {
            printf("failed to allocate memory for `P_data' (%d bytes).\n", bytes);
            endrun(1);
        }
        bytes_tot += All.MaxPart * sizeof (struct particle_data);

        P = P_data - 1; /* start with offset 1 */
        Part = &P[1];
        printf("Allocated %g MByte for particle storage.\n", bytes_tot / (1024.0 * 1024.0));
    }

}

/* This routine frees the memory for the particle storage,
 */
void free_memory(void) {

    if (All.MaxPart > 0) {
        free(P_data);
    }

}

struct classSortByEst {

    bool operator() (const int& lhs, const int& rhs) const {
        return Part[lhs].Est > Part[rhs].Est;
    }// here offset is 0 eg: index starting from 0
};

/* GADGET 2 Format for reading Ngb and Est */
void read_ic12_est(const char *fname) {
#define GETBLKNAME fread(&blkname, sizeof(blkname), 1, fd);cout<<blkname.name[0]<<blkname.name[1]<<blkname.name[2]<<blkname.name[3]<<endl;
#define SKIP my_fread1(&blklen,patterns,1,All.flag_swap,fd);
    FILE *fd;
    int i;

    int patternd[] = {4, 3, 0}, patterns[] = {4, 1, 0}, patternh[] = {4, 6, 8, 8, 4, 10, 8, 4, 4, 3, 4, 21, 0};
    int4byte blklen;
    int NumPart = 0, ONLY_TYPE = All.ONLY_TYPE; // this will select only specific type of the particles

    All.flag_swap = 0;


    if ((fd = fopen(fname, "rb"))) {
        fprintf(stdout, "Reading GADGET format file: %s \n", fname);
        fflush(stdout);


        SKIP;
        if (blklen != 256) {
            All.flag_swap = 1;
            SwapEndian(&blklen, patterns);
            if (blklen != 256) {
                printf("incorrect header format (2)\n");
                endrun(889);
            }
        }

        my_fread1(&header1, patternh, 1, All.flag_swap, fd);


        SKIP;
        if (blklen != 256) {
            printf("incorrect header format (2)\n");
            endrun(889);
        }

        All.DesNumNgb = (int) header1.BoxSize;
        All.DesNumNgbA = (int) header1.time;

        for (i = 0, NumPart = 0; i < 6; i++)
            NumPart += header1.npart[i];
        if (NumPart != All.NumPart) {
            cout << "\nError:\n\nsomething is odd Est file and Snap file are not compartable..." << endl;
            endrun(13);
        }


        SKIP;
        float *vec = new float[NumPart];
        int nread;
        nread = fread((char*) & vec[0], sizeof (float), NumPart, fd);
        if (nread != NumPart)
            endrun(778);
        for (i = 1; i <= NumPart; i++) {
            P[i].Est = vec[i - 1];
        }
        SKIP;
        delete [] vec;
        int *veci = new int[All.DesNumNgb];
        SKIP;
        for (i = 1; i <= NumPart; i++) {
            nread = fread((char*) & veci[0], sizeof (int), All.DesNumNgb, fd);
            P[i].pNGB = new int[All.DesNumNgb];
            memcpy(&P[i].pNGB[0], &veci[0], All.DesNumNgb * sizeof (int));
            std::sort(&P[i].pNGB[0], &P[i].pNGB[All.DesNumNgb - 1], classSortByEst());
            /*for(int ingb=0;ingb<All.DesNumNgb;ingb++)
                    {
                    cout<<P[i].pNGB[ingb]<<" "<<Part[P[i].pNGB[ingb]].Est<<endl;
                    }
             */
        }
        SKIP;
        delete []veci;
        if (All.DesNumNgbA > 0) {
            int *vec = new int[All.DesNumNgbA];
            SKIP;
            for (i = 1; i <= NumPart; i++) {
                nread = fread((char*) & vec[0], sizeof (int), All.DesNumNgbA, fd);
                P[i].pNGBA = new int[All.DesNumNgbA];
                memcpy(&P[i].pNGBA[0], &vec[0], All.DesNumNgbA * sizeof (int));
                std::sort(&P[i].pNGBA[0], &P[i].pNGBA[All.DesNumNgbA - 1], classSortByEst());
                /*for(int ingb=0;ingb<All.DesNumNgb;ingb++)
                {
                cout<<P[i].pNGB[ingb]<<" "<<Part[P[i].pNGB[ingb]].Est<<endl;
                }
                 */
            }
            SKIP;
            delete []vec;

        }

        fclose(fd);
        fprintf(stdout, "....done with reading.\n");
        fflush(stdout);


    } else {
        fprintf(stdout, "File %s not found.\n", fname);
        endrun(7);
    }



#undef SKIP    
#undef GETBLKNAME  
}
//End for GADGET 2 format for Est

/* GADGET 2 Format for reading Rho*/
void read_ic12_rho(const char *fname) {
#define GETBLKNAME fread(&blkname, sizeof(blkname), 1, fd);cout<<blkname.name[0]<<blkname.name[1]<<blkname.name[2]<<blkname.name[3]<<endl;
#define SKIP my_fread1(&blklen,patterns,1,All.flag_swap,fd);
    FILE *fd;
    int i;

    struct {
        int blks;
        char name[4];
        int next;
        int blk2;
    } blkname;

    int patternd[] = {4, 3, 0}, patterns[] = {4, 1, 0}, patternh[] = {4, 6, 8, 8, 4, 10, 8, 4, 4, 3, 4, 21, 0};
    int4byte blklen;
    int NumPart = 0, ONLY_TYPE = All.ONLY_TYPE; // this will select only specific type of the particles

    All.flag_swap = 0;


    if ((fd = fopen(fname, "rb"))) {
        fprintf(stdout, "Reading GADGET format file: %s \n", fname);
        fflush(stdout);

        GETBLKNAME;
        SKIP;
        if (blklen != 256) {
            All.flag_swap = 1;
            SwapEndian(&blklen, patterns);
            if (blklen != 256) {
                printf("incorrect header format (2)\n");
                endrun(889);
            }
        }

        my_fread1(&header1, patternh, 1, All.flag_swap, fd);


        SKIP;
        if (blklen != 256) {
            printf("incorrect header format (2)\n");
            endrun(889);
        }

        for (i = 0, NumPart = 0; i < 6; i++)
            NumPart += header1.npart[i];
        if (NumPart != All.NumPart) {
            cout << "\nError:\n\nsomething is odd Est file and Snap file are not compartable..." << endl;
            endrun(13);
        }

        float *loc_vec = new float[NumPart];
        int nread;

        GETBLKNAME; // ID
        SKIP;
        nread = fread((char*) & loc_vec[0], sizeof (float), NumPart, fd);
        SKIP;

        GETBLKNAME; //  RHO
        SKIP;
        nread = fread((char*) & loc_vec[0], sizeof (float), NumPart, fd);
        SKIP;
        for (i = 1; i <= NumPart; i++) {
            P[i].Rho = loc_vec[i - 1];
        }
        GETBLKNAME; // HSML
        SKIP;
        nread = fread((char*) & loc_vec[0], sizeof (float), NumPart, fd);
        if (nread != NumPart)
            endrun(778);
        for (i = 1; i <= NumPart; i++) {
            P[i].Hsml = loc_vec[i - 1];
        }
        SKIP;
        delete [] loc_vec;


        fclose(fd);
        fprintf(stdout, "....done with reading.\n");
        fflush(stdout);


    } else {
        fprintf(stdout, "File %s not found.\n", fname);
        endrun(7);
    }



#undef SKIP    
#undef GETBLKNAME  
}
//End for GADGET 2 format for Rho

double mydrand48() {
    return double(rand()) / RAND_MAX;
}

#undef NDIM
#undef real
#undef real1

