
#include "header.h"
#include "faces.h"
#include "geometry.h"
#include "hydro.h"
#include "plm_aos.h"
#include "profiler.h"
#include "substep_aos.h"


void recon_aos(struct domain *theDomain)
{
    prof_tick(theDomain->prof, PROF_RECON_P);
    plm_phi_aos(theDomain);
    prof_tock(theDomain->prof, PROF_RECON_P);

    if(theDomain->Nr > 1)
    {
        prof_tick(theDomain->prof, PROF_RECON_R);
        setup_faces( theDomain , 1 );
        int Nfr = theDomain->fIndex_r[theDomain->N_ftracks_r];
        plm_trans_aos(theDomain, theDomain->theFaces_1, Nfr, 1);
        prof_tock(theDomain->prof, PROF_RECON_R);
    }
    if(theDomain->Nz > 1)
    {
        prof_tick(theDomain->prof, PROF_RECON_Z);
        setup_faces( theDomain , 2 );
        int Nfz = theDomain->fIndex_z[theDomain->N_ftracks_z];
        plm_trans_aos(theDomain, theDomain->theFaces_2, Nfz, 2);
        prof_tock(theDomain->prof, PROF_RECON_Z);
    }
}

void addFlux_aos(struct domain *theDomain)
{
    prof_tick(theDomain->prof, PROF_FLUX_P);
    flux_phi_aos(theDomain);
    prof_tock(theDomain->prof, PROF_FLUX_P);

    if(theDomain->Nr > 1)
    {
        prof_tick(theDomain->prof, PROF_FLUX_R);
        flux_r_aos(theDomain);
        prof_tock(theDomain->prof, PROF_FLUX_R);
    }
    if(theDomain->Nz > 1)
    {
        prof_tick(theDomain->prof, PROF_FLUX_Z);
        flux_z_aos(theDomain);
        prof_tock(theDomain->prof, PROF_FLUX_Z);
    }
}

void addSource_aos(struct domain *theDomain)
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
                struct cell *c = &(theDomain->theCells[jk][i]);
                double phip = c->piph;
                double phim = c->piph - c->dphi;
                double xp[3] = {rp, phip, zp};
                double xm[3] = {rm, phim, zm};
                double dV = get_dV(xp, xm);

                source(c->prim, c->cons, xp, xm, 1.0e-3 * dV);
                visc_source(c->prim, c->gradr, c->gradp, c->gradz, c->cons,
                            xp, xm, 1.0e-3 * dV);
            }
        }
    }
}

void calcPrim_aos(struct domain *theDomain)
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
                struct cell *c = &(theDomain->theCells[jk][i]);
                double phip = c->piph;
                double phim = c->piph - c->dphi;
                double xp[3] = {rp, phip, zp};
                double xm[3] = {rm, phim, zm};
                double phi = get_centroid(phip, phim, 0);
                double dV = get_dV(xp, xm);
                double x[3] = {r, phi, z};

                cons2prim(c->cons, c->prim, x, dV, xp, xm);
            }
        }
    }
}


void flux_phi_aos(struct domain *theDomain)
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

                struct cell *cL = &(theDomain->theCells[jk][iL]);
                struct cell *cR = &(theDomain->theCells[jk][iR]);

                riemann_phi_aos(cL, cR, 1.0e-3, rm, rp, r, zm, zp, z);
            }
        }
    }
}

void flux_r_aos(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;

    int jmin, jmax, kmin, kmax, Nfr;

    int *fI;
    struct face * theFaces;

    fI = theDomain->fIndex_r;
    theFaces = theDomain->theFaces_1;
    Nfr = Nr-1;
    jmin = NgRa==0 ? 0 : NgRa-1;
    jmax = NgRb==0 ? Nr-1 : Nr-NgRb;
    kmin = NgZa;
    kmax = Nz-NgZb;

    int j, k;
    for(k=kmin; k<kmax; k++)
    {
        double zm, zp;
        zm = theDomain->z_kph[k-1];
        zp = theDomain->z_kph[k];

        for(j=jmin; j<jmax; j++)
        {
            double r = theDomain->r_jph[j];

            int JK = j + Nfr*k;
            int f;
            for(f=fI[JK]; f<fI[JK+1]; f++)
                riemann_r_aos(theFaces + f, 1.0e-3, r, zp, zm);
        }
    }
}
void flux_z_aos(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;
    
    int jmin, jmax, kmin, kmax, Nfr;

    int *fI;
    struct face * theFaces;

    fI = theDomain->fIndex_z;
    theFaces = theDomain->theFaces_2;
    Nfr = Nr;
    jmin = NgRa;
    jmax = Nr-NgRb;
    kmin = NgZa==0 ? 0  : NgZa-1;
    kmax = NgZb==0 ? Nz-1 : Nz-NgZb;


    int j, k;
    for(k=kmin; k<kmax; k++)
    {
        double z = theDomain->z_kph[k];

        for(j=jmin; j<jmax; j++)
        {
            double rm, rp;
            rm = theDomain->r_jph[j-1];
            rp = theDomain->r_jph[j];

            int JK = j + Nfr*k;
            int f;
            for(f=fI[JK]; f<fI[JK+1]; f++)
                riemann_z_aos(theFaces + f, 1.0e-3, rp, rm, z);
        }
    }
}

void riemann_phi_aos(struct cell *cL, struct cell *cR, double dt,
                     double rm, double rp, double r, double zm, double zp,
                     double z)
{   
    double xp[3] = {rp, cL->piph, zp};
    double xm[3] = {rm, cL->piph, zm};
    double x[3] = {r, cL->piph, z};
    double dA = get_dA(xp,xm,0); 

    double primL[NUM_Q];
    double primR[NUM_Q];
    double prim[NUM_Q];

    double dphiL = 0.5*cL->dphi;
    double dphiR = 0.5*cR->dphi;
    
    double n[3] = {0.0, 1.0, 0.0};

    int q;
    for(q=0; q<NUM_Q; q++)
    {
        primL[q] = cL->prim[q] + dphiL * cL->gradp[q];
        primR[q] = cR->prim[q] - dphiR * cR->gradp[q];
        prim[q] = 0.5*(primL[q] + primR[q]);
    }

    double FL[NUM_Q], FR[NUM_Q], F[NUM_Q], VF[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        VF[q] = 0.0;

    flux(primL, FL, x, n, xp, xm);
    flux(primR, FR, x, n, xp, xm);
    visc_flux(prim, cL->gradr, cL->gradp, cL->gradz,
              VF, x, n);
    
    for(q=0; q<NUM_Q; q++)
        F[q] = 0.5*(FL[q] + FR[q]);

    for(q=0; q<NUM_Q; q++)
    {
        cL->cons[q] -= (F[q] + VF[q]) * dt * dA;
        cR->cons[q] += (F[q] + VF[q]) * dt * dA;
    }
}
                
void riemann_r_aos(struct face *f, double dt, double r, double zp, double zm)
{
    double primL[NUM_Q], primR[NUM_Q], prim[NUM_Q];

    double n[3] = {1.0, 0.0, 0.0};
    double xp[3] = {r, f->cm[1] + 0.5*f->dphi, zp};
    double xm[3] = {r, f->cm[1] - 0.5*f->dphi, zm};

    double dphiL = f->cm[1] - (f->L->piph - 0.5*f->L->dphi);
    double dphiR = f->cm[1] - (f->R->piph - 0.5*f->R->dphi);

    int q;
    for(q=0; q<NUM_Q; q++)
    {
        primL[q] = f->L->prim[q] + dphiL * f->L->gradp[q]
                    + f->dxL * f->L->gradr[q];
        primR[q] = f->R->prim[q] + dphiR * f->R->gradp[q]
                    + f->dxR * f->R->gradr[q];
        prim[q] = 0.5*(primL[q] + primR[q]);
    }

    double FL[NUM_Q], FR[NUM_Q], F[NUM_Q], VF[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        VF[q] = 0.0;

    flux(primL, FL, f->cm, n, xp, xm);
    flux(primR, FR, f->cm, n, xp, xm);
    visc_flux(prim, f->L->gradr, f->L->gradp, f->L->gradz,
              VF, f->cm, n);
    
    for(q=0; q<NUM_Q; q++)
        F[q] = 0.5*(FL[q] + FR[q]);

    for(q=0; q<NUM_Q; q++)
    {
        f->L->cons[q] -= (F[q] + VF[q]) * dt * f->dA * dt;
        f->R->cons[q] += (F[q] + VF[q]) * dt * f->dA * dt;
    }
}
                
void riemann_z_aos(struct face *f, double dt, double rp, double rm, double z)
{
    double primL[NUM_Q], primR[NUM_Q], prim[NUM_Q];

    double n[3] = {0.0, 0.0, 1.0};
    double xp[3] = {rp, f->cm[1] + 0.5*f->dphi, z};
    double xm[3] = {rm, f->cm[1] - 0.5*f->dphi, z};

    double dphiL = f->cm[1] - (f->L->piph - 0.5*f->L->dphi);
    double dphiR = f->cm[1] - (f->R->piph - 0.5*f->R->dphi);

    int q;
    for(q=0; q<NUM_Q; q++)
    {
        primL[q] = f->L->prim[q] + dphiL * f->L->gradp[q]
                    + f->dxL * f->L->gradz[q];
        primR[q] = f->R->prim[q] + dphiR * f->R->gradp[q]
                    + f->dxR * f->R->gradz[q];
        prim[q] = 0.5*(primL[q] + primR[q]);
    }

    double FL[NUM_Q], FR[NUM_Q], F[NUM_Q], VF[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        VF[q] = 0.0;

    flux(primL, FL, f->cm, n, xp, xm);
    flux(primR, FR, f->cm, n, xp, xm);
    visc_flux(prim, f->L->gradr, f->L->gradp, f->L->gradz,
              VF, f->cm, n);
    
    for(q=0; q<NUM_Q; q++)
        F[q] = 0.5*(FL[q] + FR[q]);
    
    for(q=0; q<NUM_Q; q++)
    {
        f->L->cons[q] -= (F[q] + VF[q]) * dt * f->dA;
        f->R->cons[q] += (F[q] + VF[q]) * dt * f->dA;
    }
}
