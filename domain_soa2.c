
#include "header.h"
#include "geometry.h"
#include "hydro.h"
#include "domain_soa2.h"

void setupDomain_soa2(struct domain *theDomain)
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
    
    int j, k;
    theDomain->fr_iL = (int **)malloc(Nr*Nz*sizeof(int *));
    theDomain->fr_dA_L = (double **)malloc(Nr*Nz*sizeof(double *));
    theDomain->fr_dphi_L = (double **)malloc(Nr*Nz*sizeof(double *));
    theDomain->fr_offset_L = (double **)malloc(Nr*Nz*sizeof(double *));
    for(k=0; k<Nz; k++)
    {
        int jk = Nr*k + 0;
        theDomain->fr_iL[jk] = NULL;
        theDomain->fr_dA_L[jk] = NULL;
        theDomain->fr_dphi_L[jk] = NULL;
        theDomain->fr_offset_L[jk] = NULL;

        for(j=1; j < Nr; j++)
        {
            jk = Nr*k + j;
            theDomain->fr_iL[jk] = (int *)malloc(Np[jk] * sizeof(int));
            theDomain->fr_dA_L[jk] = (double *)malloc(Np[jk] * sizeof(double));
            theDomain->fr_dphi_L[jk] = (double *)malloc(Np[jk]
                                                        * sizeof(double));
            theDomain->fr_offset_L[jk] = (double *)malloc(Np[jk] * sizeof(double));
        }
    }
    theDomain->fr_iR = (int **)malloc(Nr*Nz*sizeof(int *));
    theDomain->fr_dA_R = (double **)malloc(Nr*Nz*sizeof(double *));
    theDomain->fr_dphi_R = (double **)malloc(Nr*Nz*sizeof(double *));
    theDomain->fr_offset_R = (double **)malloc(Nr*Nz*sizeof(double *));
    for(k=0; k<Nz; k++)
    {
        int jk = Nr*k + Nr-1;
        theDomain->fr_iR[jk] = NULL;
        theDomain->fr_dA_R[jk] = NULL;
        theDomain->fr_dphi_R[jk] = NULL;
        theDomain->fr_offset_R[jk] = NULL;

        for(j=0; j < Nr-1; j++)
        {
            jk = Nr*k + j;
            theDomain->fr_iR[jk] = (int *)malloc(Np[jk] * sizeof(int));
            theDomain->fr_dA_R[jk] = (double *)malloc(Np[jk] * sizeof(double));
            theDomain->fr_dphi_R[jk] = (double *)malloc(Np[jk]
                                                        * sizeof(double));
            theDomain->fr_offset_R[jk] = (double *)malloc(Np[jk] * sizeof(double));
        }
    }

    theDomain->fz_iL = (int **)malloc(Nr*Nz*sizeof(int *));
    theDomain->fz_dA_L = (double **)malloc(Nr*Nz*sizeof(double *));
    theDomain->fz_dphi_L = (double **)malloc(Nr*Nz*sizeof(double *));
    theDomain->fz_offset_L = (double **)malloc(Nr*Nz*sizeof(double *));
    for(j=0; j < Nr; j++)
    {
        int jk = Nr*(0) + j;
        theDomain->fz_iL[jk] = NULL;
        theDomain->fz_dA_L[jk] = NULL;
        theDomain->fz_dphi_L[jk] = NULL;
        theDomain->fz_offset_L[jk] = NULL;
    }

    for(k=1; k<Nz; k++)
        for(j=0; j < Nr; j++)
        {
            int jk = Nr*k + j;
            theDomain->fz_iL[jk] = (int *)malloc(Np[jk] * sizeof(int));
            theDomain->fz_dA_L[jk] = (double *)malloc(Np[jk] * sizeof(double));
            theDomain->fz_dphi_L[jk] = (double *)malloc(Np[jk]
                                                        * sizeof(double));
            theDomain->fz_offset_L[jk] = (double *)malloc(Np[jk] * sizeof(double));
        }

    theDomain->fz_iR = (int **)malloc(Nr*Nz*sizeof(int *));
    theDomain->fz_dA_R = (double **)malloc(Nr*Nz*sizeof(double *));
    theDomain->fz_dphi_R = (double **)malloc(Nr*Nz*sizeof(double *));
    theDomain->fz_offset_R = (double **)malloc(Nr*Nz*sizeof(double *));
    
    for(j=0; j < Nr; j++)
    {
        int jk = Nr*(Nz-1) + j;
        theDomain->fz_iR[jk] = NULL;
        theDomain->fz_dA_R[jk] = NULL;
        theDomain->fz_dphi_R[jk] = NULL;
        theDomain->fz_offset_R[jk] = NULL;
    }

    for(k=0; k<Nz; k++)
        for(j=0; j < Nr; j++)
        {
            int jk = Nr*k + j;
            theDomain->fz_iR[jk] = (int *)malloc(Np[jk] * sizeof(int));
            theDomain->fz_dA_R[jk] = (double *)malloc(Np[jk] * sizeof(double));
            theDomain->fz_dphi_R[jk] = (double *)malloc(Np[jk]
                                                        * sizeof(double));
            theDomain->fz_offset_R[jk] = (double *)malloc(Np[jk] * sizeof(double));
        }

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


    if(theDomain->Nr > 1)
        build_faces_soa2(theDomain, DIM_R);
    if(theDomain->Nz > 1)
        build_faces_soa2(theDomain, DIM_Z);

    check_faces(theDomain);
}

void build_faces_soa2(struct domain *theDomain, int dim)
{
    if(dim != DIM_R && dim != DIM_Z)
        return;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    double phi_max = theDomain->phi_max;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;
    
    int k;
    for(k=0; k<Nz; k++)
    {
        int j;
        for(j=0; j<Nr; j++)
        {
            int jk = Nr*k + j;
            int jkL, jkR;

            if(dim == DIM_R)
            {
                jkL = jk - 1;
                jkR = jk + 1;
                
                if(j > 0 && theDomain->fr_iL[jk] != NULL)
                {
                    fill_faces_soa2(Np[jk], Np[jkL], piph[jk], piph[jkL],
                         dphi[jk], dphi[jkL],
                         r_jph[j-1], r_jph[j-1], z_kph[k-1], z_kph[k],
                         theDomain->fr_iL[jk], theDomain->fr_dphi_L[jk],
                         theDomain->fr_dA_L[jk], theDomain->fr_offset_L[jk],
                         phi_max, 1, dim);
                }

                if(j < Nr-1 && theDomain->fr_iR[jk] != NULL)
                {
                    fill_faces_soa2(Np[jk], Np[jkR], piph[jk], piph[jkR],
                         dphi[jk], dphi[jkR],
                         r_jph[j], r_jph[j], z_kph[k-1], z_kph[k],
                         theDomain->fr_iR[jk], theDomain->fr_dphi_R[jk],
                         theDomain->fr_dA_R[jk], theDomain->fr_offset_R[jk],
                         phi_max, 1, dim);
                }

            }
            else
            {
                jkL = jk - Nr;
                jkR = jk + Nr;

                if(k > 0 && theDomain->fz_iL[jk] != NULL)
                {
                    fill_faces_soa2(Np[jk], Np[jkL], piph[jk], piph[jkL],
                         dphi[jk], dphi[jkL],
                         r_jph[j-1], r_jph[j], z_kph[k-1], z_kph[k-1],
                         theDomain->fz_iL[jk], theDomain->fz_dphi_L[jk],
                         theDomain->fz_dA_L[jk], theDomain->fz_offset_L[jk],
                         phi_max, 1, dim);
                }
                if(k < Nz-1 && theDomain->fz_iR[jk] != NULL)
                {
                    fill_faces_soa2(Np[jk], Np[jkR], piph[jk], piph[jkR],
                         dphi[jk], dphi[jkR],
                         r_jph[j-1], r_jph[j], z_kph[k], z_kph[k],
                         theDomain->fz_iR[jk], theDomain->fz_dphi_R[jk],
                         theDomain->fz_dA_R[jk], theDomain->fz_offset_R[jk],
                         phi_max, 1, dim);
                }
            }
        }
    }
}

void fill_faces_soa2(int Np, int Np1, const double *piph, const double *piph1,
                     const double *dphi, const double *dphi1,
                     double rm, double rp, double zm, double zp,
                     int *f_i, double *f_dphi, double *f_dA, double *f_offset,
                     double phi_max, int inc_match, int dim)
{
    int i1 = 0;
    int i;

    while(get_signed_dp(piph1[i1], piph[0]) < 0.0 
            || get_signed_dp(piph1[i1] - dphi1[i1], piph[0]) > 0.0)
    {
        i1 = (i1 == Np1-1) ? 0 : i1+1;
    }

    for(i=0; i<Np; i++)
    {
        // compute the offset between cells i and i1;
        double offset = 0.0;
        while((piph1[i1]+offset) - piph[i] > 0.5*phi_max)
            offset -= phi_max;
        while((piph1[i1]+offset) - piph[i] < -0.5*phi_max)
            offset += phi_max;

        // front of cell i1 adjusted to be on i's branch
        double phif1 = piph1[i1] + offset;

        // iterate i1 forwards until it intesects i's front
        while(phif1 < piph[i] || phif1 - dphi1[i1] > piph[i])
        {
            i1 = (i1 == Np1-1) ? 0 : i1+1;

            // update offset
            while((piph1[i1]+offset) - piph[i] > 0.5*phi_max)
                offset -= phi_max;
            while((piph1[i1]+offset) - piph[i] < -0.5*phi_max)
                offset += phi_max;

            phif1 = piph1[i1] + offset;
        }

        // TODO: handle match case for quad

        //Find back of face
        double phib = phif1 - dphi1[i1];
        if(phib < piph[i] - dphi[i])
            phib = piph[i] - dphi[i];
        double phif = piph[i];

        double xm[3] = {rm, phib, zm};
        double xp[3] = {rp, phif, zp};

        f_i[i] = i1;
        f_dphi[i] = phif-phib;
        f_dA[i] = get_dA(xp, xm, dim);
        f_offset[i] = phif1 - phif;
    }
}

void check_faces(struct domain *theDomain)
{
    int Nz = theDomain->Nz;
    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;

    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int **fr_iL = theDomain->fr_iL;
    int **fr_iR = theDomain->fr_iR;

    int k, j, i;
    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Nr; j++)
        {
            int jk = k*Nr + j;
            int jkL = jk-1;

            if(fr_iL[jk] == NULL)
                continue;

            for(i=0; i<Np[jk]; i++)
            {
                int iL = fr_iL[jk][i];
                double phiL = piph[jkL][iL];
                double dphiL = dphi[jkL][iL];

                double phiR = piph[jk][i];

                if(get_signed_dp(phiL, phiR) < 0.0)
                {
                    printf("BAD FACE rL too back\n");
                    continue;
                }
                if(get_signed_dp(phiL-dphiL, phiR) > 0.0)
                {
                    printf("BAD FACE rL too fore\n");
                    continue;
                }

                double phibR = phiR - dphi[jk][i];
                double phibL = phiL - dphi[jkL][iL];

                double phif = get_signed_dp(phiL, phiR) > 0.0 ? phiR : phiL;
                double phib = get_signed_dp(phibL, phibR) > 0.0 ? phibL : phibR;
                double r = theDomain->r_jph[j-1];
                double zm = theDomain->z_kph[k-1];
                double zp = theDomain->z_kph[k];

                double xp[3] = {r, phif, zp};
                double xm[3] = {r, phib, zm};

                double dA = get_dA(xp, xm, 1);
                double dphif = get_signed_dp(phif, phib);
                double offset = get_signed_dp(phiL, phiR);

                double dAf = theDomain->fr_dA_L[jk][i];
                double dphiff = theDomain->fr_dphi_L[jk][i];
                double offsetf = theDomain->fr_offset_L[jk][i];

                if(fabs(dA - dAf) > 1.0e-14)
                    printf("BAD FACE wrong dA  %lg %lg %le\n",
                            dA, dAf, dA-dAf);
                if(fabs(dphif - dphiff) > 1.0e-14)
                    printf("BAD FACE wrong phib  %lg %lg %le\n",
                            dphif, dphiff, dphif-dphiff);
                if(fabs(offset - offsetf) > 1.0e-14)
                    printf("BAD FACE wrong offset  %lg %lg %le\n",
                            offset, offsetf, offset-offsetf);

            }
        }
    }
    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Nr; j++)
        {
            int jk = k*Nr + j;
            int jkR = jk+1;

            if(fr_iR[jk] == NULL)
                continue;

            for(i=0; i<Np[jk]; i++)
            {
                int iR = fr_iR[jk][i];
                double phiR = piph[jkR][iR];
                double dphiR = dphi[jkR][iR];

                double phiL = piph[jk][i];

                if(get_signed_dp(phiR, phiL) < 0.0)
                {
                    printf("BAD FACE rR too back\n");
                    continue;
                }
                if(get_signed_dp(phiR-dphiR, phiL) > 0.0)
                {
                    printf("BAD FACE rR too fore\n");
                    continue;
                }

                double phibL = phiL - dphi[jk][i];
                double phibR = phiR - dphi[jkR][iR];

                double phif = get_signed_dp(phiL, phiR) > 0.0 ? phiR : phiL;
                double phib = get_signed_dp(phibL, phibR) > 0.0 ? phibL : phibR;
                double r = theDomain->r_jph[j];
                double zm = theDomain->z_kph[k-1];
                double zp = theDomain->z_kph[k];

                double xp[3] = {r, phif, zp};
                double xm[3] = {r, phib, zm};

                double dA = get_dA(xp, xm, 1);
                double dphif = get_signed_dp(phif, phib);
                double offset = get_signed_dp(phiR, phiL);

                double dAf = theDomain->fr_dA_R[jk][i];
                double dphiff = theDomain->fr_dphi_R[jk][i];
                double offsetf = theDomain->fr_offset_R[jk][i];

                if(fabs(dA - dAf) > 1.0e-14)
                    printf("BAD FACE wrong dA  %lg %lg %le\n",
                            dA, dAf, dA-dAf);
                if(fabs(dphif - dphiff) > 1.0e-14)
                    printf("BAD FACE wrong dphif  %lg %lg %le\n",
                            dphif, dphiff, dphif-dphiff);
                if(fabs(offset - offsetf) > 1.0e-14)
                    printf("BAD FACE wrong offset  %lg %lg %le\n",
                            offset, offsetf, offset-offsetf);
            }
        }
    }
}

void freeDomain_soa2(struct domain *theDomain)
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
        if(theDomain->fr_iL[jk] != NULL)
            free(theDomain->fr_iL[jk]);
        if(theDomain->fr_dA_L[jk] != NULL)
            free(theDomain->fr_dA_L[jk]);
        if(theDomain->fr_dphi_L[jk] != NULL)
            free(theDomain->fr_dphi_L[jk]);
        if(theDomain->fr_offset_L[jk] != NULL)
            free(theDomain->fr_offset_L[jk]);
        if(theDomain->fr_iR[jk] != NULL)
            free(theDomain->fr_iR[jk]);
        if(theDomain->fr_dA_R[jk] != NULL)
            free(theDomain->fr_dA_R[jk]);
        if(theDomain->fr_dphi_R[jk] != NULL)
            free(theDomain->fr_dphi_R[jk]);
        if(theDomain->fr_offset_R[jk] != NULL)
            free(theDomain->fr_offset_R[jk]);
        if(theDomain->fz_iL[jk] != NULL)
            free(theDomain->fz_iL[jk]);
        if(theDomain->fz_dA_L[jk] != NULL)
            free(theDomain->fz_dA_L[jk]);
        if(theDomain->fz_dphi_L[jk] != NULL)
            free(theDomain->fz_dphi_L[jk]);
        if(theDomain->fz_offset_L[jk] != NULL)
            free(theDomain->fz_offset_L[jk]);
        if(theDomain->fz_iR[jk] != NULL)
            free(theDomain->fz_iR[jk]);
        if(theDomain->fz_dA_R[jk] != NULL)
            free(theDomain->fz_dA_R[jk]);
        if(theDomain->fz_dphi_R[jk] != NULL)
            free(theDomain->fz_dphi_R[jk]);
        if(theDomain->fz_offset_R[jk] != NULL)
            free(theDomain->fz_offset_R[jk]);
    }
    free(theDomain->prim);
    free(theDomain->cons);
    free(theDomain->gradr);
    free(theDomain->gradp);
    free(theDomain->gradz);
    free(theDomain->piph);
    free(theDomain->dphi);

    free(theDomain->fr_iL);
    free(theDomain->fr_dA_L);
    free(theDomain->fr_dphi_L);
    free(theDomain->fr_offset_L);
    free(theDomain->fr_iR);
    free(theDomain->fr_dA_R);
    free(theDomain->fr_dphi_R);
    free(theDomain->fr_offset_R);
    free(theDomain->fz_iL);
    free(theDomain->fz_dA_L);
    free(theDomain->fz_dphi_L);
    free(theDomain->fz_offset_L);
    free(theDomain->fz_iR);
    free(theDomain->fz_dA_R);
    free(theDomain->fz_dphi_R);
    free(theDomain->fz_offset_R);

    free(theDomain->Np);
    theDomain->r_jph--;
    free(theDomain->r_jph);
    theDomain->z_kph--;
    free(theDomain->z_kph);
}

double hash_soa2(struct domain *theDomain, int qqq)
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

void dump_grid_soa2(struct domain *theDomain, char label[])
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
