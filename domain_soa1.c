
#include "header.h"
#include "geometry.h"
#include "hydro.h"
#include "domain_soa1.h"

void setupDomain_soa1(struct domain *theDomain)
{
    srand(314159);
    rand();

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;

    int jk;

    // Allocate everything
    theDomain->prim = (double **)malloc(Nr*Nz*sizeof(double *));
    for(jk=0; jk < Nr*Nz; jk++)
        theDomain->prim[jk] = (double *)malloc(Np[jk]*NUM_Q*sizeof(double *));
    
    theDomain->cons = (double **)malloc(Nr*Nz*sizeof(double *));
    for(jk=0; jk < Nr*Nz; jk++)
        theDomain->cons[jk] = (double *)malloc(Np[jk]*NUM_Q*sizeof(double *));
    
    theDomain->gradr = (double **)malloc(Nr*Nz*sizeof(double *));
    for(jk=0; jk < Nr*Nz; jk++)
        theDomain->gradr[jk] = (double *)malloc(Np[jk]*NUM_Q*sizeof(double *));
    
    theDomain->gradp = (double **)malloc(Nr*Nz*sizeof(double *));
    for(jk=0; jk < Nr*Nz; jk++)
        theDomain->gradp[jk] = (double *)malloc(Np[jk]*NUM_Q*sizeof(double *));
    
    theDomain->gradz = (double **)malloc(Nr*Nz*sizeof(double *));
    for(jk=0; jk < Nr*Nz; jk++)
        theDomain->gradz[jk] = (double *)malloc(Np[jk]*NUM_Q*sizeof(double *));
    
    theDomain->dphi = (double **)malloc(Nr*Nz*sizeof(double *));
    for(jk=0; jk < Nr*Nz; jk++)
        theDomain->dphi[jk] = (double *)malloc(Np[jk]*sizeof(double *));
    
    theDomain->piph = (double **)malloc(Nr*Nz*sizeof(double *));
    for(jk=0; jk < Nr*Nz; jk++)
        theDomain->piph[jk] = (double *)malloc(Np[jk]*sizeof(double *));


    // Set phi face locations
    int i, j, k;
    double Pmax = theDomain->phi_max;

    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = Nr*k + j;
            double p0 = Pmax*(double)rand()/(double)RAND_MAX;
            double dp = Pmax/(double)Np[jk];
            for( i=0 ; i<Np[jk] ; ++i )
            {
                double phi = p0 + dp*i;
                if( phi > Pmax )
                    phi -= Pmax;
                theDomain->piph[jk][i] = phi;
                theDomain->dphi[jk][i] = dp;
            }

        }

    // Initialize hydro data;
    for(k=0; k<Nz; k++)
    {
        double zp = theDomain->z_kph[k];
        double zm = theDomain->z_kph[k-1];
        double z = get_centroid(zp, zm, 2);

        for(j=0; j<Nr; j++)
        {
            double rp = theDomain->r_jph[j];
            double rm = theDomain->r_jph[j-1];
            double r = get_centroid(rp, rm, 1);
            
            int jk = Nr*k + j;
            for(i=0; i<Np[jk]; i++)
            {
                double phip = theDomain->piph[jk][i];
                double phim = phip - theDomain->dphi[jk][i];
                double phi = 0.5*(phim + phip);

                double xp[3] = {rp, phip, zp};
                double xm[3] = {rm, phim, zm};
                double x[3] = {r, phi, z};

                double dV = get_dV(xp, xm);

                int iq = i * NUM_Q;

                initial(theDomain->prim[jk]+iq, x);
                prim2cons(theDomain->prim[jk]+iq, theDomain->cons[jk]+iq,
                          x, dV, xp, xm);
                cons2prim(theDomain->cons[jk]+iq, theDomain->prim[jk]+iq,
                          x, dV, xp, xm);

                int q;
                for(q=0; q<NUM_Q; q++)
                {
                    theDomain->gradr[jk][iq+q] = 0.0;
                    theDomain->gradp[jk][iq+q] = 0.0;
                    theDomain->gradz[jk][iq+q] = 0.0;
                }
            }
        }
    }
}


void freeDomain_soa1(struct domain *theDomain)
{
    int jk;
    int Nann = theDomain->Nr * theDomain->Nz;

    for(jk=0; jk<Nann; jk++)
    {
        free(theDomain->prim[jk]);
        free(theDomain->cons[jk]);
        free(theDomain->gradr[jk]);
        free(theDomain->gradp[jk]);
        free(theDomain->gradz[jk]);
        free(theDomain->piph[jk]);
        free(theDomain->dphi[jk]);
    }
    free(theDomain->prim);
    free(theDomain->cons);
    free(theDomain->gradr);
    free(theDomain->gradp);
    free(theDomain->gradz);
    free(theDomain->piph);
    free(theDomain->dphi);
}

double hash_soa1(struct domain *theDomain, int qqq)
{
    int i, j, k;
    double sum = 0;

    for(k=0; k<theDomain->Nz; k++)
        for(j=0; j<theDomain->Nr; j++)
        {
            int jk = k * theDomain->Nr + j;
            for(i=0; i<theDomain->Np[jk]; i++)
            {
                if(qqq >= 0 && qqq < NUM_Q)
                    sum += theDomain->prim[jk][i*NUM_Q+qqq]; 
                else
                {
                    int q;
                    for(q=0; q<NUM_Q; q++)
                        sum += theDomain->prim[jk][i*NUM_Q+q]; 
                }
            }
        }

    return sum;
}

