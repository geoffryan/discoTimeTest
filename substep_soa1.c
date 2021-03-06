#include "header.h"
#include "geometry.h"
#include "hydro.h"
#include "plm_soa1.h"
#include "profiler.h"
#include "substep_soa1.h"

void recon_soa1(struct domain *theDomain)
{
    prof_tick(theDomain->prof, PROF_RECON_P);
    plm_phi_soa1(theDomain);
    prof_tock(theDomain->prof, PROF_RECON_P);

    /*

    if(theDomain->Nr > 1)
    {
        prof_tick(theDomain->prof, PROF_RECON_R);
        setup_faces( theDomain , 1 );
        int Nfr = theDomain->fIndex_r[theDomain->N_ftracks_r];
        plm_trans_soa1(theDomain, theDomain->theFaces_1, Nfr, 1);
        prof_tock(theDomain->prof, PROF_RECON_R);
    }
    if(theDomain->Nz > 1)
    {
        prof_tick(theDomain->prof, PROF_RECON_Z);
        setup_faces( theDomain , 2 );
        int Nfz = theDomain->fIndex_z[theDomain->N_ftracks_z];
        plm_trans_soa1(theDomain, theDomain->theFaces_2, Nfz, 2);
        prof_tock(theDomain->prof, PROF_RECON_Z);
    }
    */
}

void addFlux_soa1(struct domain *theDomain)
{
    prof_tick(theDomain->prof, PROF_FLUX_P);
    flux_phi_soa1(theDomain);
    prof_tock(theDomain->prof, PROF_FLUX_P);
}

void addSource_soa1(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;

    int k;
    for(k=NgZa; k<Nz-NgZb; k++)
    {
        double zm = theDomain->z_kph[k-1];
        double zp = theDomain->z_kph[k];
        int j;
        for(j=NgRa; j<Nr-NgRb; j++)
        {
            double rm = theDomain->r_jph[j-1];
            double rp = theDomain->r_jph[j];

            int jk = k*Nr + j;
            int i;
            for(i=0; i<Np[jk]; i++)
            {
                double phip = theDomain->piph[jk][i];
                double phim = phip - theDomain->dphi[jk][i];
                double xp[3] = {rp, phip, zp};
                double xm[3] = {rm, phim, zm};
                double dV = get_dV(xp, xm);

                int iq = NUM_Q*i;
                double *prim = theDomain->prim[jk] + iq;
                double *cons = theDomain->cons[jk] + iq;
                double *gradr = theDomain->gradr[jk] + iq;
                double *gradp = theDomain->gradp[jk] + iq;
                double *gradz = theDomain->gradz[jk] + iq;

                source(prim, cons, xp, xm, 1.0e-3 * dV);
                visc_source(prim, gradr, gradp, gradz, cons,
                            xp, xm, 1.0e-3 * dV);
            }
        }
    }
}

void calcPrim_soa1(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;

    int k;
    for(k=NgZa; k<Nz-NgZb; k++)
    {
        double zm = theDomain->z_kph[k-1];
        double zp = theDomain->z_kph[k];
        double z = get_centroid(zp, zm, 2);
        int j;
        for(j=NgRa; j<Nr-NgRb; j++)
        {
            double rm = theDomain->r_jph[j-1];
            double rp = theDomain->r_jph[j];
            double r = get_centroid(rp, rm, 1);

            int jk = k*Nr + j;
            int i;
            for(i=0; i<Np[jk]; i++)
            {
                double phip = theDomain->piph[jk][i];
                double phim = phip - theDomain->dphi[jk][i];
                double xp[3] = {rp, phip, zp};
                double xm[3] = {rm, phim, zm};
                double x[3] = {r, 0.5*(phim+phip), z};

                double dV = get_dV(xp, xm);

                int iq = NUM_Q*i;
                double *prim = theDomain->prim[jk] + iq;
                double *cons = theDomain->cons[jk] + iq;

                cons2prim(cons, prim, x, dV, xp, xm);
            }
        }
    }
}

void flux_phi_soa1(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;

    int k;
    for(k=NgZa; k<Nz-NgZb; k++)
    {
        double zm = theDomain->z_kph[k-1];
        double zp = theDomain->z_kph[k];
        double z = get_centroid(zp, zm, 2);
        int j;
        for(j=NgRa; j<Nr-NgRb; j++)
        {
            double rm = theDomain->r_jph[j-1];
            double rp = theDomain->r_jph[j];
            double r = get_centroid(rp, rm, 1);

            int jk = k*Nr + j;
            int iL;
            for(iL=0; iL<Np[jk]; iL++)
            {
                int iR = iL < Np[jk]-1 ? iL + 1 : 0;

                riemann_phi_soa1(theDomain, jk, iL, iR, 1.0e-3,
                                 rm, rp, r, zm, zp, z);
            }
        }
    }
    
}

void riemann_phi_soa1(struct domain *theDomain, int jk, int iL, int iR,
                      double dt, double rm, double rp, double r,
                      double zm, double zp, double z)
{   
    double xp[3] = {rp, theDomain->piph[jk][iL], zp};
    double xm[3] = {rm, theDomain->piph[jk][iL], zm};
    double x[3] = {r, theDomain->piph[jk][iL], z};
    double dA = get_dA(xp,xm,0); 

    double primL[NUM_Q];
    double primR[NUM_Q];
    double prim[NUM_Q];

    double dphiL = 0.5*theDomain->dphi[jk][iL];
    double dphiR = 0.5*theDomain->dphi[jk][iR];
    
    double n[3] = {0.0, 1.0, 0.0};

    int iqL = NUM_Q * iL;
    int iqR = NUM_Q * iR;

    int q;
    for(q=0; q<NUM_Q; q++)
    {
        primL[q] = theDomain->prim[jk][iqL+q]
                    + dphiL * theDomain->gradp[jk][iqL+q];
        primR[q] = theDomain->prim[jk][iqR+q]
                    - dphiR * theDomain->gradp[jk][iqR+q];
    }
    for(q=0; q<NUM_Q; q++)
        prim[q] = 0.5*(primL[q] + primR[q]);

    double F[NUM_Q], VF[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        VF[q] = 0.0;

    flux(prim, F, x, n, xp, xm);
    visc_flux(prim, theDomain->gradr[jk]+iqL, theDomain->gradp[jk]+iqL,
              theDomain->gradz[jk]+iqL, VF, x, n);

    for(q=0; q<NUM_Q; q++)
    {
        theDomain->cons[jk][iqL+q] -= (F[q] + VF[q]) * dt * dA;
        theDomain->cons[jk][iqR+q] += (F[q] + VF[q]) * dt * dA;
    }
}
