#include "header.h"
#include "domain_soa2.h"
#include "geometry.h"
#include "hydro.h"
#include "plm_soa2.h"
#include "profiler.h"
#include "substep_soa2.h"

void recon_soa2(struct domain *theDomain)
{
    prof_tick(theDomain->prof, PROF_RECON_P);
    plm_phi_soa2(theDomain);
    prof_tock(theDomain->prof, PROF_RECON_P);


    if(theDomain->Nr > 1)
    {
        prof_tick(theDomain->prof, PROF_RECON_R);
        build_faces_soa2(theDomain);
        plm_r_soa2(theDomain);
        prof_tock(theDomain->prof, PROF_RECON_R);
    }
    /*
    if(theDomain->Nz > 1)
    {
        prof_tick(theDomain->prof, PROF_RECON_Z);
        build_faces_soa2(theDomain);
        plm_z_soa2(theDomain);
        prof_tock(theDomain->prof, PROF_RECON_Z);
    }
    */
}

void addFlux_soa2(struct domain *theDomain, double dt)
{
    prof_tick(theDomain->prof, PROF_FLUX_P);
    flux_phi_soa2(theDomain, dt);
    prof_tock(theDomain->prof, PROF_FLUX_P);
    prof_tick(theDomain->prof, PROF_FLUX_R);
    flux_r_soa2(theDomain, dt);
    prof_tock(theDomain->prof, PROF_FLUX_R);
}

void addSource_soa2(struct domain *theDomain, double dt)
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

                source(prim, cons, xp, xm, dt * dV);
                visc_source(prim, gradr, gradp, gradz, cons,
                            xp, xm, dt * dV);
            }
        }
    }
}

void calcPrim_soa2(struct domain *theDomain)
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
                double x[3] = {r, get_centroid(phip, phim, 0), z};

                double dV = get_dV(xp, xm);

                int iq = NUM_Q*i;
                double *prim = theDomain->prim[jk] + iq;
                double *cons = theDomain->cons[jk] + iq;

                cons2prim(cons, prim, x, dV, xp, xm);
            }
        }
    }
}

void flux_phi_soa2(struct domain *theDomain, double dt)
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

#if SUBTYPE == 0
                riemann_phi_soa2(theDomain, jk, iL, iR, dt,
                                 rm, rp, r, zm, zp, z);
#elif SUBTYPE == 1
                riemann_phi_soa2_alt(theDomain->prim[jk], theDomain->cons[jk],
                                     theDomain->dphi[jk], theDomain->piph[jk],
                                     theDomain->gradr[jk], theDomain->gradp[jk],
                                     theDomain->gradz[jk],
                                     iL, iR, dt,
                                     rm, rp, r, zm, zp, z);
#endif
            }
        }
    }
    
}

void riemann_phi_soa2(struct domain *theDomain, int jk, int iL, int iR,
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
        prim[q] = 0.5*(primL[q] + primR[q]);
    }

    double FL[NUM_Q], FR[NUM_Q], F[NUM_Q], VF[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        VF[q] = 0.0;

    flux(primL, FL, x, n, xp, xm);
    flux(primR, FR, x, n, xp, xm);
    visc_flux(prim, theDomain->gradr[jk]+iqL, theDomain->gradp[jk]+iqL,
              theDomain->gradz[jk]+iqL, VF, x, n);
    for(q=0; q<NUM_Q; q++)
        F[q] = 0.5*(FL[q] + FR[q]);

    for(q=0; q<NUM_Q; q++)
    {
        theDomain->cons[jk][iqL+q] -= (F[q] + VF[q]) * dt * dA;
        theDomain->cons[jk][iqR+q] += (F[q] + VF[q]) * dt * dA;
    }
}

void riemann_phi_soa2_alt(const double *prim, double *cons,
                          const double *dphi, const double *piph,
                          const double *gradr, const double *gradp,
                          const double *gradz,
                          int iL, int iR,
                          double dt, double rm, double rp, double r,
                          double zm, double zp, double z)
{   
    double xp[3] = {rp, piph[iL], zp};
    double xm[3] = {rm, piph[iL], zm};
    double x[3] = {r, piph[iL], z};
    double dA = get_dA(xp,xm,0); 

    double primL[NUM_Q];
    double primR[NUM_Q];
    double primF[NUM_Q];

    double dphiL = 0.5*dphi[iL];
    double dphiR = 0.5*dphi[iR];
    
    double n[3] = {0.0, 1.0, 0.0};

    int iqL = NUM_Q * iL;
    int iqR = NUM_Q * iR;

    int q;
    for(q=0; q<NUM_Q; q++)
    {
        primL[q] = prim[iqL+q] + dphiL * gradp[iqL+q];
        primR[q] = prim[iqR+q] - dphiR * gradp[iqR+q];
    }
    for(q=0; q<NUM_Q; q++)
        primF[q] = 0.5*(primL[q] + primR[q]);

    double F[NUM_Q], VF[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        VF[q] = 0.0;

    flux(primF, F, x, n, xp, xm);
    visc_flux(primF, gradr+iqL, gradp+iqL, gradz+iqL, VF, x, n);

    for(q=0; q<NUM_Q; q++)
    {
        cons[iqL+q] -= (F[q] + VF[q]) * dt * dA;
        cons[iqR+q] += (F[q] + VF[q]) * dt * dA;
    }
}

void flux_r_soa2(struct domain *theDomain, double dt)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;

    double *r_jph = theDomain->r_jph;
    int *Nfr = theDomain->Nfr;

    int k;
    for(k=NgZa; k<Nz-NgZb; k++)
    {
        double zm = theDomain->z_kph[k-1];
        double zp = theDomain->z_kph[k];
        double z = get_centroid(zp, zm, 2);
        int j;

        for(j=NgRa-1; j<Nr-NgRb; j++)
        {
            double rL = get_centroid(r_jph[j], r_jph[j-1], 1);
            double rR = get_centroid(r_jph[j+1], r_jph[j], 1);
            double r = r_jph[j];

            int jkf = (Nr-1)*k + j;
            int jkL = k*Nr + j;
            int jkR = k*Nr + j+1;

            double *piphL = theDomain->piph[jkL];
            double *piphR = theDomain->piph[jkR];

#if SUBTYPE == 1
            double *primL = theDomain->prim[jkL];
            double *consL = theDomain->cons[jkL];
            double *gradrL = theDomain->gradr[jkL];
            double *gradpL = theDomain->gradp[jkL];
            double *gradzL = theDomain->gradz[jkL];
            double *dphiL = theDomain->dphi[jkL];
            
            double *primR = theDomain->prim[jkR];
            double *consR = theDomain->cons[jkR];
            double *gradrR = theDomain->gradr[jkR];
            double *gradpR = theDomain->gradp[jkR];
            double *gradzR = theDomain->gradz[jkR];
            double *dphiR = theDomain->dphi[jkR];
            
            double *fr_dA = theDomain->fr_dA[jkf];
            double *fr_phif = theDomain->fr_phif[jkf];
            double *fr_phib = theDomain->fr_phib[jkf];
#endif

            int i;
            int iL = theDomain->I0[jkL];
            int iR = theDomain->I0[jkR];
            for(i=0; i<Nfr[jkf]; i++)
            {
#if SUBTYPE == 0
                riemann_r_soa2(theDomain, jkL, jkR, jkf, iL, iR, i,
                                 dt, rL, rR, r, zm, zp, z);
#elif SUBTYPE == 1
                riemann_r_soa2_alt(primL, primR, consL, consR,
                                   gradrL, gradpL, gradzL,
                                   gradrR, gradpR, gradzR,
                                   dphiL, piphL, dphiR, piphR,
                                   fr_dA, fr_phif, fr_phib,
                                   iL, iR, i, dt, rL, rR, r, zm, zp, z);
#endif

                double dpLR = get_signed_dp(piphL[iL], piphR[iR]);
                if(dpLR < 0)
                    iL = iL==Np[jkL]-1 ? 0 : iL+1;
                else
                    iR = iR==Np[jkR]-1 ? 0 : iR+1;
            }
        }
    }
}

void riemann_r_soa2(struct domain *theDomain, int jkL, int jkR, int jkf,
                    int iL, int iR, int i, double dt, double rL, double rR,
                    double r, double zm, double zp, double z)
{
    double xp[3] = {r, theDomain->fr_phif[jkf][i], zp};
    double xm[3] = {r, theDomain->fr_phib[jkf][i], zm};
    double dA = theDomain->fr_dA[jkf][i];

    double primL[NUM_Q];
    double primR[NUM_Q];
    double prim[NUM_Q];

    double phiF = 0.5*(theDomain->fr_phif[jkf][i]
                       + theDomain->fr_phib[jkf][i]);
    double phiL = theDomain->piph[jkL][iL] - 0.5*theDomain->dphi[jkL][iL];
    double phiR = theDomain->piph[jkR][iR] - 0.5*theDomain->dphi[jkR][iR];

    double dphiL = get_signed_dp(phiF, phiL);
    double dphiR = get_signed_dp(phiF, phiR);
    double drL = r - rL;
    double drR = r - rR;
    double x[3] = {r, phiF, z};
    
    double n[3] = {1.0, 0.0, 0.0};

    int iqL = NUM_Q * iL;
    int iqR = NUM_Q * iR;

    int q;
    for(q=0; q<NUM_Q; q++)
    {
        primL[q] = theDomain->prim[jkL][iqL+q]
                    + dphiL * theDomain->gradp[jkL][iqL+q]
                    + drL * theDomain->gradr[jkL][iqL+q];
        primR[q] = theDomain->prim[jkR][iqR+q]
                    + dphiR * theDomain->gradp[jkR][iqR+q]
                    + drR * theDomain->gradr[jkR][iqR+q];
        prim[q] = 0.5*(primL[q] + primR[q]);
    }

    double FL[NUM_Q], FR[NUM_Q], F[NUM_Q], VF[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        VF[q] = 0.0;

    flux(primL, FL, x, n, xp, xm);
    flux(primR, FR, x, n, xp, xm);
    for(q=0; q<NUM_Q; q++)
        F[q] = 0.5*(FL[q] + FR[q]);
    
    visc_flux(prim, theDomain->gradr[jkL]+iqL, theDomain->gradp[jkL]+iqL,
              theDomain->gradz[jkL]+iqL, VF, x, n);

    for(q=0; q<NUM_Q; q++)
    {
        theDomain->cons[jkL][iqL+q] -= (F[q] + VF[q]) * dt * dA;
        theDomain->cons[jkR][iqR+q] += (F[q] + VF[q]) * dt * dA;
    }
}

void riemann_r_soa2_alt(const double *primL, const double *primR,
                        double *consL, double *consR,
                        const double *gradrL, const double *gradpL,
                        const double *gradzL,
                        const double *gradrR, const double *gradpR,
                        const double *gradzR,
                        const double *dphiL, const double *piphL,
                        const double *dphiR, const double *piphR,
                        const double *fr_dA, const double *fr_phif,
                        const double *fr_phib,
                        int iL, int iR, int i, double dt, double rL, double rR,
                        double r, double zm, double zp, double z)
{
    double xp[3] = {r, fr_phif[i], zp};
    double xm[3] = {r, fr_phib[i], zm};
    double dA = fr_dA[i];

    double pL[NUM_Q];
    double pR[NUM_Q];
    double prim[NUM_Q];

    double phiF = 0.5*(xp[1] + xm[1]);
    double phiL = piphL[iL] - 0.5*dphiL[iL];
    double phiR = piphR[iR] - 0.5*dphiR[iR];

    double dpL = get_signed_dp(phiF, phiL);
    double dpR = get_signed_dp(phiF, phiR);
    double drL = r - rL;
    double drR = r - rR;
    double x[3] = {r, phiF, z};
    
    double n[3] = {1.0, 0.0, 0.0};

    int iqL = NUM_Q * iL;
    int iqR = NUM_Q * iR;

    int q;
    for(q=0; q<NUM_Q; q++)
    {
        pL[q] = primL[iqL+q] + dpL * gradpL[iqL+q] + drL * gradrL[iqL+q];
        pR[q] = primR[iqR+q] + dpR * gradpR[iqR+q] + drR * gradrR[iqR+q];
        prim[q] = 0.5*(pL[q] + pR[q]);
    }

    double FL[NUM_Q], FR[NUM_Q], F[NUM_Q], VF[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        VF[q] = 0.0;

    flux(pL, FL, x, n, xp, xm);
    flux(pR, FR, x, n, xp, xm);
    for(q=0; q<NUM_Q; q++)
        F[q] = 0.5*(FL[q] + FR[q]);
    
    visc_flux(prim, gradrL+iqL, gradpL+iqL, gradzL+iqL, VF, x, n);

    for(q=0; q<NUM_Q; q++)
    {
        consL[iqL+q] -= (F[q] + VF[q]) * dt * dA;
        consR[iqR+q] += (F[q] + VF[q]) * dt * dA;
    }
}
