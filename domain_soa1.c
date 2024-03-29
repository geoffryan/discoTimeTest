
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
    
    theDomain->I0 = (int *)malloc(Nr*Nz*sizeof(int));
    
    int j, k;
    theDomain->Nfr = (int *)malloc((Nr-1)*Nz*sizeof(int));
    for(k=0; k<Nz; k++)
        for(j=0; j < Nr-1; j++)
            theDomain->Nfr[(Nr-1)*k+j] = Np[k*Nr+j] + Np[k*Nr+j+1];
    int *Nfr = theDomain->Nfr;
    theDomain->fr_dA = (double **)malloc((Nr-1)*Nz * sizeof(double *));
    for(jk=0; jk < (Nr-1)*Nz; jk++)
        theDomain->fr_dA[jk] = (double *)malloc(Nfr[jk] * sizeof(double));
    theDomain->fr_phib = (double **)malloc((Nr-1)*Nz * sizeof(double *));
    for(jk=0; jk < (Nr-1)*Nz; jk++)
        theDomain->fr_phib[jk] = (double *)malloc(Nfr[jk] * sizeof(double));
    theDomain->fr_phif = (double **)malloc((Nr-1)*Nz * sizeof(double *));
    for(jk=0; jk < (Nr-1)*Nz; jk++)
        theDomain->fr_phif[jk] = (double *)malloc(Nfr[jk] * sizeof(double));
    
    theDomain->Nfz = (int *)malloc(Nr*(Nz-1)*sizeof(int));
    for(k=0; k<Nz-1; k++)
        for(j=0; j < Nr; j++)
            theDomain->Nfz[Nr*k+j] = Np[k*Nr+j] + Np[(k+1)*Nr+j];
    int *Nfz = theDomain->Nfz;
    theDomain->fz_dA = (double **)malloc(Nr*(Nz-1) * sizeof(double *));
    for(jk=0; jk < Nr*(Nz-1); jk++)
        theDomain->fz_dA[jk] = (double *)malloc(Nfz[jk] * sizeof(double));
    theDomain->fz_phib = (double **)malloc(Nr*(Nz-1) * sizeof(double *));
    for(jk=0; jk < Nr*(Nz-1); jk++)
        theDomain->fz_phib[jk] = (double *)malloc(Nfz[jk] * sizeof(double));
    theDomain->fz_phif = (double **)malloc(Nr*(Nz-1) * sizeof(double *));
    for(jk=0; jk < Nr*(Nz-1); jk++)
        theDomain->fz_phif[jk] = (double *)malloc(Nfz[jk] * sizeof(double));
    

    // Set phi face locations
    int i;
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

    //reset dphi
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = Nr*k + j;
            for( i=0 ; i<Np[jk] ; ++i )
            {
                int im = i == 0 ? Np[jk] - 1 : i-1;
                double phip = theDomain->piph[jk][i];
                double phim = theDomain->piph[jk][im];
                theDomain->dphi[jk][i] = get_dp(phip, phim);
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


    build_faces_soa1(theDomain);
}

void build_faces_soa1(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    double phi_max = theDomain->phi_max;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;
    int *I0 = theDomain->I0;

    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *Nfr = theDomain->Nfr;
    double **fr_dA = theDomain->fr_dA;
    double **fr_phib = theDomain->fr_phib;
    double **fr_phif = theDomain->fr_phif;
    
    int *Nfz = theDomain->Nfz;
    double **fz_dA = theDomain->fz_dA;
    double **fz_phib = theDomain->fz_phib;
    double **fz_phif = theDomain->fz_phif;

    int jk;
    
    // Find faces crossing phi=0 & re-align piph to (0, phi_max]
    for(jk=0; jk<Nr*Nz; jk++)
    {
        int i;

        while(piph[jk][Np[jk]-1] > phi_max)
            piph[jk][Np[jk]-1] -= phi_max;
        while(piph[jk][Np[jk]-1] <= 0)
            piph[jk][Np[jk]-1] += phi_max;
        
        I0[jk] = -1;

        int adjust_check = 0;

        for(i=0; i<Np[jk]; i++)
        {
            while(piph[jk][i] > phi_max)
                piph[jk][i] -= phi_max;
            while(piph[jk][i] <= 0)
                piph[jk][i] += phi_max;

            int im = i == 0 ? Np[jk]-1 : i-1;
            dphi[jk][i] = piph[jk][i] - piph[jk][im];

            if(dphi[jk][i] < 0.0)
            {
                I0[jk] = i;
                dphi[jk][i] += phi_max;
                adjust_check++;
            }
        }

        if(I0[jk] == -1)
            printf("BAD ALIGNMENT\n");
        if(adjust_check != 1)
            printf("Bad # of adjustments: %d\n", adjust_check);
    }

    // Now that we have a set of adjacent cells, we can find & compute faces

    // radial faces
    int k;
    for(k=0; k<Nz; k++)
    {
        double zm = z_kph[k-1];
        double zp = z_kph[k];
        
        int j;
        for(j=0; j<Nr-1; j++)
        {
            int jkL = k*Nr + j;
            int jkR = k*Nr + j+1;
            int jkf = k*(Nr-1) + j;

            double r = r_jph[j];

            int iL = I0[jkL];
            int iR = I0[jkR];
            int i;
            for(i=0; i<Nfr[jkf]; i++)
            {
                double phifL = piph[jkL][iL];
                double phifR = piph[jkR][iR];
                if(phifR - phifL > 0.5*phi_max)
                    phifL += phi_max;
                else if(phifR - phifL < -0.5*phi_max)
                    phifR += phi_max;

                double phibL = phifL - dphi[jkL][iL];
                double phibR = phifR - dphi[jkR][iR];
                double phif = phifL > phifR ? phifR : phifL;
                double phib = phibL > phibR ? phibL : phibR;
                
                //printf("    %d %d %d   %.03lf %.03lf\n", i, iL, iR,
                //        phifL, phifR);

                double xp[3] = {r, phif, zp};
                double xm[3] = {r, phib, zm};

                fr_dA[jkf][i] = get_dA(xp, xm, 1);
                fr_phib[jkf][i] = phib;
                fr_phif[jkf][i] = phif;

                double dpLR = get_signed_dp(phifL, phifR);

                if(dpLR < 0.0)
                    iL = iL == Np[jkL]-1 ? 0 : iL+1;
                else
                    iR = iR == Np[jkR]-1 ? 0 : iR+1;
            }
            if(iL != I0[jkL] || iR != I0[jkR])
                printf("Radial faces didn't finish: %d %d (%d) %d %d (%d)\n",
                        iL, I0[jkL], Np[jkL], iR, I0[jkR], Np[jkR]);
        }
    }

    for(k=0; k<Nz-1; k++)
    {
        double z = z_kph[k];
        
        int j;
        for(j=0; j<Nr; j++)
        {
            int jkL = k*Nr + j;
            int jkR = (k+1)*Nr + j;
            int jkf = k*Nr + j;

            double rm = r_jph[j-1];
            double rp = r_jph[j];

            int iL = I0[jkL];
            int iR = I0[jkR];
            int i;
            for(i=0; i<Nfz[jkf]; i++)
            {
                double phifL = piph[jkL][iL];
                double phifR = piph[jkR][iR];
                double phibL = phifL - dphi[jkL][iL];
                double phibR = phifR - dphi[jkR][iR];
                double phif = phifL > phifR ? phifR : phifL;
                double phib = phibL > phibR ? phibL : phibR;

                double xp[3] = {rp, phif, z};
                double xm[3] = {rm, phib, z};

                fz_dA[jkf][i] = get_dA(xp, xm, 2);
                fz_phib[jkf][i] = phib;
                fz_phif[jkf][i] = phif;

                double dpLR = phifL - phifR;
                while(dpLR > 0.5*phi_max)
                    dpLR -= phi_max;
                while(dpLR < -0.5*phi_max)
                    dpLR += phi_max;

                if(dpLR < 0.0)
                {
                    iL++;
                    if(iL == Np[jkL])
                        iL = 0;
                }
                else
                {
                    iR++;
                    if(iR == Np[jkR])
                        iR = 0;
                }
            }
            if(iL != I0[jkL] || iR != I0[jkR])
                printf("Vertical faces didn't finish: %d %d (%d) %d %d (%d)\n",
                        iL, I0[jkL], Np[jkL], iR, I0[jkR], Np[jkR]);
        }
    }
}


void freeDomain_soa1(struct domain *theDomain)
{
    int jk;
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int Nann = Nr * Nz;

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

    free(theDomain->Np);
    theDomain->r_jph--;
    free(theDomain->r_jph);
    theDomain->z_kph--;
    free(theDomain->z_kph);

    free(theDomain->I0);

    for(jk=0; jk<(Nr-1)*Nz; jk++)
    {
        free(theDomain->fr_dA[jk]);
        free(theDomain->fr_phib[jk]);
        free(theDomain->fr_phif[jk]);
    }
    free(theDomain->fr_dA);
    free(theDomain->fr_phib);
    free(theDomain->fr_phif);

    for(jk=0; jk<Nr*(Nz-1); jk++)
    {
        free(theDomain->fz_dA[jk]);
        free(theDomain->fz_phib[jk]);
        free(theDomain->fz_phif[jk]);
    }
    free(theDomain->fz_dA);
    free(theDomain->fz_phib);
    free(theDomain->fz_phif);

    free(theDomain->Nfr);
    free(theDomain->Nfz);
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
                    {
                        sum += theDomain->prim[jk][i*NUM_Q+q]; 
                        //sum += theDomain->gradr[jk][i*NUM_Q+q]; 
                        //sum += theDomain->gradp[jk][i*NUM_Q+q]; 
                        //sum += theDomain->gradz[jk][i*NUM_Q+q]; 
                    }
                }
            }
        }

    return sum;
}

void dump_grid_soa1(struct domain *theDomain, char label[])
{
    int i, j, k;

    char filename[512];
    sprintf(filename, "grid.%d%d.%03d.%s.txt", TYPE, SUBTYPE,
            theDomain->current_step, label);
    FILE *f = fopen(filename, "w");

    for(k=0; k<theDomain->Nz; k++)
        for(j=0; j<theDomain->Nr; j++)
        {
            int jk = j + theDomain->Nr * k;

            for(i=0; i<theDomain->Np[jk]; i++)
            {
                fprintf(f, "%03d %03d %04d:", k, j, i);

                int q;
                for(q = 0; q < NUM_Q; q++)
                    fprintf(f, " %.12le", theDomain->prim[jk][i*NUM_Q+q]);
                fprintf(f, "   %.12le", theDomain->piph[jk][i]);
                fprintf(f, "   %.12le", theDomain->dphi[jk][i]);
                fprintf(f, "\n              ");
                for(q = 0; q < NUM_Q; q++)
                    fprintf(f, " %.12le", theDomain->cons[jk][i*NUM_Q+q]);
                fprintf(f, "\n              ");
                for(q = 0; q < NUM_Q; q++)
                    fprintf(f, " %.12le", theDomain->gradr[jk][i*NUM_Q+q]);
                fprintf(f, "\n              ");
                for(q = 0; q < NUM_Q; q++)
                    fprintf(f, " %.12le", theDomain->gradp[jk][i*NUM_Q+q]);
                fprintf(f, "\n              ");
                for(q = 0; q < NUM_Q; q++)
                    fprintf(f, " %.12le", theDomain->gradz[jk][i*NUM_Q+q]);
                fprintf(f, "\n");
            }
        }
}
