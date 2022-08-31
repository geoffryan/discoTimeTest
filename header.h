#ifndef HEADER_H
#define HEADER_H

enum{RHO,PPP,URR,UPP,UZZ,BRR,BPP,BZZ};
enum{DDD,TAU,SRR,LLL,SZZ};
enum{PLPOINTMASS, PLPW, PLSURFACEGRAV, PLALEX, PLQUAD};
enum{DIM_P, DIM_R, DIM_Z};

enum{PROF_TOT, PROF_DT, PROF_TIMESTEP, PROF_OUTPUT, PROF_RECON, PROF_FLUX,
     PROF_CT, PROF_SOURCE, PROF_C2P, PROF_BOUND, PROF_EXCHANGE,
     PROF_RECON_R, PROF_RECON_P, PROF_RECON_Z, 
     PROF_FLUX_R, PROF_FLUX_P, PROF_FLUX_Z,
     PROF_EXCH_NP_COUNT1, PROF_EXCH_NP_COMM1, PROF_EXCH_NP_COUNT2,
     PROF_EXCH_NP_COMM2, PROF_EXCH_NP_FIN,
     PROF_EXCH_PREP, PROF_EXCH_COMM, PROF_EXCH_FIN,
     NUM_PROF}; // NUM_PROF must be at end

#if USE_MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// NUM_C, NUM_N, and CT_MODE are specified at compile time and defined
// with -D in the Makefile

#define NUM_Q (NUM_C+NUM_N)
#define NUM_G 2

#if CT_MODE == 0        //No CT
    #define NUM_EDGES 0    
    #define NUM_FACES 0    
    #define NUM_AZ_EDGES 0 
#elif CT_MODE == 1      //3D MHD
    #define NUM_EDGES 8    
    #define NUM_FACES 5    
    #define NUM_AZ_EDGES 4 
#endif

struct param_list{

   int nsteps;
   int Num_R, Num_Z;

   double rmin, rmax;
   double zmin, zmax;
   double phimax;
   double dt;
};

struct domain{ 
   struct profiler *prof;
   int * Np;
   int Nr,Nz,Ng;
   int NgRa, NgRb, NgZa, NgZb;
   int N0r, N0z;
   int N_ftracks_r;
   int N_ftracks_z;
   int Npl;
   double * r_jph;
   double * z_kph;
   double phi_max;

   time_t Wallt_init;

   struct param_list theParList;

   int current_step;
   int final_step;

   // AOS setup
   struct cell ** theCells;
   struct face * theFaces_1;
   struct face * theFaces_2;
   int * fIndex_r;
   int * fIndex_z;

   //SOA1 setup
   double **prim;
   double **cons;
   double **gradr;
   double **gradp;
   double **gradz;
   double **dphi;
   double **piph;

   int *I0;

   int *Nfr;
   double **fr_dA;
   double **fr_phib;
   double **fr_phif;

   int *Nfz;
   double **fz_dA;
   double **fz_phib;
   double **fz_phif;

   //SOA2 setup
   int **fr_iL;
   double **fr_dA_L;
   double **fr_phib_L;
   double **fr_phif_L;
   double **fr_dphiL;

   int **fr_iR;
   double **fr_dA_R;
   double **fr_phib_R;
   double **fr_phif_R;
   double **fr_dphiR;
   
   int **fz_iL;
   double **fz_dA_L;
   double **fz_phib_L;
   double **fz_phif_L;
   double **fz_dphiL;

   int **fz_iR;
   double **fz_dA_R;
   double **fz_phib_R;
   double **fz_phif_R;
   double **fz_dphiR;
};

struct cell{

   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double gradr[NUM_Q];
   double gradp[NUM_Q];
   double gradz[NUM_Q];
   double piph;
   double dphi;
   double wiph;

   double E[NUM_EDGES];
   double B[NUM_EDGES];
   double E_phi[NUM_AZ_EDGES];
   double    Phi[NUM_FACES];
   double RK_Phi[NUM_FACES];
   double tempDoub;

   int real;
};

struct edge{
   struct cell * LU;
   struct cell * RU;
   struct cell * LD;
   struct cell * RD;

   int Prim_Zone;
   int Alt_LR;
   int Alt_UD;

   double E_dl;
};

struct face{
   struct cell * L;
   struct cell * R;
   double dxL;
   double dxR;
   double cm[3];
   double dphi;
   double dl;
   double dA;

   double E,B;
   int LRtype;
   int flip_flag;
};

struct profiler{
    clock_t ticks[NUM_PROF];
    clock_t elapsed_ticks[NUM_PROF];
};

#endif
