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
        build_faces_soa2(theDomain, DIM_R);
        plm_r_soa2(theDomain);
        prof_tock(theDomain->prof, PROF_RECON_R);
    }
    /*
    if(theDomain->Nz > 1)
    {
        prof_tick(theDomain->prof, PROF_RECON_Z);
        build_faces_soa2(theDomain, DIM_Z);
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

                riemann_phi_soa2(theDomain, jk, iL, iR, dt,
                                 rm, rp, r, zm, zp, z);
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

void flux_r_soa2(struct domain *theDomain, double dt)
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
            
            const double *piph = theDomain->piph[jk];
            const double *dphi = theDomain->dphi[jk];

            const double *prim = theDomain->prim[jk];
            const double *gradr = theDomain->gradr[jk];
            const double *gradp = theDomain->gradp[jk];
            const double *gradz = theDomain->gradz[jk];
            double *cons = theDomain->cons[jk];

            if(j > 0)
            {
                int jkL = jk - 1;
                const double *dphiL = theDomain->dphi[jkL];
                double rmm = theDomain->r_jph[j-2];
                double rL = get_centroid(rm, rmm, 1);
                
                const double *primL = theDomain->prim[jkL];
                const double *gradrL = theDomain->gradr[jkL];
                const double *gradpL = theDomain->gradp[jkL];
                const double *gradzL = theDomain->gradz[jkL];
                double *consL = theDomain->cons[jkL];

                const int *f_iL = theDomain->fr_iL[jk];
                const double *f_dA_L = theDomain->fr_dA_L[jk];
                const double *f_phib_L = theDomain->fr_phib_L[jk];
                const double *f_dphiL = theDomain->fr_dphiL[jk];

                int i;
                for(i=0; i<Np[jk]; i++)
                {
                    int iL = f_iL[i];
                    double x[3] = {r, piph[i]-0.5*dphi[i], z};
                    double xL[3] = {rL, piph[i]+f_dphiL[i]-0.5*dphiL[iL], z};
                    double xm[3] = {rm, f_phib_L[i], zm};
                    double xp[3] = {rm, piph[i], zp};

                    int iqL = iL * NUM_Q;
                    int iq = i * NUM_Q;

                    riemann_soa2(primL+iqL, prim+iq,
                                    consL+iqL, cons+iq,
                                    gradrL+iqL, gradpL+iqL, gradzL+iqL,
                                    gradr+iq, gradp+iq, gradz+iq,
                                    xL, x, dt, f_dA_L[i], xp, xm, DIM_R);
                }

            }

            if(j < Nr-1)
            {
                int jkR = jk + 1;
                const double *dphiR = theDomain->dphi[jkR];
                double rpp = theDomain->r_jph[j+1];
                double rR = get_centroid(rp, rpp, 1);
                
                const double *primR = theDomain->prim[jkR];
                const double *gradrR = theDomain->gradr[jkR];
                const double *gradpR = theDomain->gradp[jkR];
                const double *gradzR = theDomain->gradz[jkR];
                double *consR = theDomain->cons[jkR];

                const int *f_iR = theDomain->fr_iR[jk];
                const double *f_dA_R = theDomain->fr_dA_R[jk];
                const double *f_phib_R = theDomain->fr_phib_R[jk];
                const double *f_dphiR = theDomain->fr_dphiR[jk];

                int i;
                for(i=0; i<Np[jk]; i++)
                {
                    int iR = f_iR[i];
                    double x[3] = {r, piph[i]-0.5*dphi[i], z};
                    double xR[3] = {rR, piph[i]+f_dphiR[i]-0.5*dphiR[iR], z};
                    double xm[3] = {rp, f_phib_R[i], zm};
                    double xp[3] = {rp, piph[i], zp};

                    int iqR = iR * NUM_Q;
                    int iq = i * NUM_Q;

                    riemann_soa2(prim+iq, primR+iqR,
                                    cons+iq, consR+iqR,
                                    gradr+iq, gradp+iq, gradz+iq,
                                    gradrR+iqR, gradpR+iqR, gradzR+iqR,
                                    x, xR, dt, f_dA_R[i], xp, xm, DIM_R);
                }
            }
        }
    }
}


void flux_z_soa2(struct domain *theDomain, double dt)
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
            
            const double *piph = theDomain->piph[jk];
            const double *dphi = theDomain->dphi[jk];

            const double *prim = theDomain->prim[jk];
            const double *gradr = theDomain->gradr[jk];
            const double *gradp = theDomain->gradp[jk];
            const double *gradz = theDomain->gradz[jk];
            double *cons = theDomain->cons[jk];

            if(k > 0)
            {
                int jkL = jk - Nr;
                const double *dphiL = theDomain->dphi[jkL];
                double zmm = theDomain->z_kph[k-2];
                double zL = get_centroid(zm, zmm, 2);
                
                const double *primL = theDomain->prim[jkL];
                const double *gradrL = theDomain->gradr[jkL];
                const double *gradpL = theDomain->gradp[jkL];
                const double *gradzL = theDomain->gradz[jkL];
                double *consL = theDomain->cons[jkL];

                const int *f_iL = theDomain->fz_iL[jk];
                const double *f_dA_L = theDomain->fz_dA_L[jk];
                const double *f_phib_L = theDomain->fz_phib_L[jk];
                const double *f_dphiL = theDomain->fz_dphiL[jk];

                int i;
                for(i=0; i<Np[jk]; i++)
                {
                    int iL = f_iL[i];
                    double x[3] = {r, piph[i]-0.5*dphi[i], z};
                    double xL[3] = {r, piph[i]+f_dphiL[i]-0.5*dphiL[iL], zL};
                    double xm[3] = {rm, f_phib_L[i], zm};
                    double xp[3] = {rp, piph[i], zm};

                    int iqL = iL * NUM_Q;
                    int iq = i * NUM_Q;

                    riemann_soa2(primL+iqL, prim+iq,
                                    consL+iqL, cons+iq,
                                    gradrL+iqL, gradpL+iqL, gradzL+iqL,
                                    gradr+iq, gradp+iq, gradz+iq,
                                    xL, x, dt, f_dA_L[i], xp, xm, DIM_Z);
                }

            }

            if(k < Nz-1)
            {
                int jkR = jk + Nr;
                const double *dphiR = theDomain->dphi[jkR];
                double zpp = theDomain->z_kph[k+1];
                double zR = get_centroid(zp, zpp, 2);
                
                const double *primR = theDomain->prim[jkR];
                const double *gradrR = theDomain->gradr[jkR];
                const double *gradpR = theDomain->gradp[jkR];
                const double *gradzR = theDomain->gradz[jkR];
                double *consR = theDomain->cons[jkR];

                const int *f_iR = theDomain->fz_iR[jk];
                const double *f_dA_R = theDomain->fz_dA_R[jk];
                const double *f_phib_R = theDomain->fz_phib_R[jk];
                const double *f_dphiR = theDomain->fz_dphiR[jk];

                int i;
                for(i=0; i<Np[jk]; i++)
                {
                    int iR = f_iR[i];
                    double x[3] = {r, piph[i]-0.5*dphi[i], z};
                    double xR[3] = {r, piph[i]+f_dphiR[i]-0.5*dphiR[iR], zR};
                    double xm[3] = {rm, f_phib_R[i], zp};
                    double xp[3] = {rp, piph[i], zp};

                    int iqR = iR * NUM_Q;
                    int iq = i * NUM_Q;

                    riemann_soa2(prim+iq, primR+iqR,
                                    cons+iq, consR+iqR,
                                    gradr+iq, gradp+iq, gradz+iq,
                                    gradrR+iqR, gradpR+iqR, gradzR+iqR,
                                    x, xR, dt, f_dA_R[i], xp, xm, DIM_Z);
                }
            }
        }
    }
}

void riemann_soa2(const double *primL, const double *primR,
                  double *consL, double *consR,
                  const double *gradrL, const double *gradpL,
                  const double *gradzL, const double *gradrR,
                  const double *gradpR, const double *gradzR,
                  const double *xL, const double *xR,
                  double dt, double dA, const double *xp, const double *xm,
                  int dim)
{
    double dxL, dxR;
    double n[3] = {0.0, 0.0, 0.0};
    if(dim == DIM_R)
    {
        dxL = 0.5*(xm[0]+xp[0]) - xL[0];
        dxR = 0.5*(xm[0]+xp[0]) - xR[0];
        n[0] = 1.0;
    }
    else
    {
        dxL = 2.5*(xm[2]+xp[2]) - xL[2];
        dxR = 2.5*(xm[2]+xp[2]) - xR[2];
        n[2] = 1.0;
    }

    const double *gradL = (dim == DIM_R) ? gradrL : gradzL;
    const double *gradR = (dim == DIM_R) ? gradrR : gradzR;

    double x[3];
    get_centroid_arr(xp, xm, x);

    double dphiL = x[1] - xL[1];
    double dphiR = x[1] - xR[1];

    double pL[NUM_Q];
    double pR[NUM_Q];
    double pA[NUM_Q];

    int q;
    for(q=0; q<NUM_Q; q++)
    {
        pL[q] = primL[q] + dphiL * gradpL[q] + dxL * gradL[q];
        pR[q] = primR[q] + dphiR * gradpR[q] + dxR * gradR[q];
        pA[q] = 0.5*(pL[q] + pR[q]);
    }

    double FL[NUM_Q], FR[NUM_Q], F[NUM_Q], VF[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        VF[q] = 0.0;

    flux(primL, FL, x, n, xp, xm);
    flux(primR, FR, x, n, xp, xm);
    for(q=0; q<NUM_Q; q++)
        F[q] = 0.5*(FL[q] + FR[q]);
    
    visc_flux(pA, gradrL, gradpL, gradzL, VF, x, n);

    for(q=0; q<NUM_Q; q++)
    {
        consL[q] -= (F[q] + VF[q]) * dt * dA;
        consR[q] += (F[q] + VF[q]) * dt * dA;
    }
}

