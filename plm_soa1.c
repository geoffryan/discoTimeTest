
#include "header.h"
#include "geometry.h"
#include "plm_soa1.h"

double minmod_soa1(double a, double b, double c)
{
    if( a*b <= 0.0  || b*c <= 0.0)
        return 0.0;

    double m = a;
    
    if( fabs(b) < fabs(m) )
        m = b;
    
    if( fabs(c) < fabs(m) )
        m = c;
   
    return m;
}

void plm_phi_soa1(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;
    double PLM = 1.5;
    int i,j,k,q;

    for( k=0 ; k<Nz ; ++k )
    {
        for( j=0 ; j<Nr ; ++j )
        {
            int jk = j+Nr*k;

            double *dphi = theDomain->dphi[jk];
            double *prim = theDomain->prim[jk];
            double *gradp = theDomain->gradp[jk];

            for( i=0 ; i<Np[jk] ; ++i )
            {
                int iL = i == 0 ? Np[jk]-1 : i-1;
                int iR = i == Np[jk]-1 ? 0 : i+1;
                
                int iqL = NUM_Q*iL;
                int iqC = NUM_Q*i;
                int iqR = NUM_Q*iR;

                for( q=0 ; q<NUM_Q ; ++q )
                {
                    double pL = prim[iqL+q];
                    double pC = prim[iqC+q];
                    double pR = prim[iqR+q];
                    double sL = pC - pL;
                    sL /= .5*( dphi[i] + dphi[iL] );
                    double sR = pR - pC;
                    sR /= .5*( dphi[iR] + dphi[i] );
                    double sC = pR - pL;
                    sC /= .5*( dphi[iL] + dphi[iR] ) + dphi[i];
                    gradp[iqC+q] = minmod_soa1( PLM*sL , sC , PLM*sR );
                }
            }
        }
    }
}

void plm_r_soa1(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *Nfr = theDomain->Nfr;

    double PLM = 1.5;
   
    // Clear radial gradients
    int i, j, k, q;
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = Nr*k + j;
            double *gradr = theDomain->gradr[jk];
            for(i=0; i<Np[jk]; i++)
                for(q=0; q<NUM_Q; q++)
                    gradr[NUM_Q*i+q] = 0.0;
        }

    // Add face-centered gradients
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr-1; j++)
        {
            int jkL = j+Nr*k;
            int jkR = j+1+Nr*k;
            int jkf = j + (Nr-1)*k;

            double *piphL = theDomain->piph[jkL];
            double *dphiL = theDomain->dphi[jkL];
            double *primL = theDomain->prim[jkL];
            double *gradrL = theDomain->gradr[jkL];
            double *gradpL = theDomain->gradp[jkL];
            
            double *piphR = theDomain->piph[jkR];
            double *dphiR = theDomain->dphi[jkR];
            double *primR = theDomain->prim[jkR];
            double *gradrR = theDomain->gradr[jkR];
            double *gradpR = theDomain->gradp[jkR];
            
            double *fr_dA = theDomain->fr_dA[jkf];
            double *fr_phif = theDomain->fr_phif[jkf];
            double *fr_phib = theDomain->fr_phib[jkf];

            double rL = get_centroid(r_jph[j], r_jph[j-1], 1);
            double rR = get_centroid(r_jph[j+1], r_jph[j], 1);
            double idr = 1.0/(rR - rL);

            int iL = theDomain->I0[jkL];
            int iR = theDomain->I0[jkR];

            for( i=0 ; i<Nfr[jkf] ; ++i )
            {
                double dpL = get_signed_dp(0.5*(fr_phif[i]-fr_phib[i]), 
                                           piphL[iL] - 0.5*dphiL[iL]);
                double dpR = get_signed_dp(0.5*(fr_phif[i]-fr_phib[i]), 
                                           piphR[iR] - 0.5*dphiR[iR]);

                for( q=0 ; q<NUM_Q ; ++q )
                {
                    double fL = primL[NUM_Q*iL+q] + dpL * gradpL[NUM_Q*iL+q];
                    double fR = primR[NUM_Q*iR+q] + dpR * gradpR[NUM_Q*iR+q];
                    double gr = (fR - fL) * fr_dA[i] * idr;
                    gradrL[NUM_Q*iL+q] += gr;
                    gradrR[NUM_Q*iR+q] += gr;
                }

                double dpLR = get_signed_dp(piphL[iL], piphR[iR]);
                if(dpLR < 0)
                    iL = iL==Np[jkL]-1 ? 0 : iL+1;
                else
                    iR = iR==Np[jkR]-1 ? 0 : iR+1;
            }
        }
   
    // Divide out by total area
    for(k=0; k<Nz; k++)
    {
        double zm = z_kph[k-1];
        double zp = z_kph[k];

        for(j=0; j<Nr; j++)
        {
            int jk = Nr*k + j;

            double rm = r_jph[j-1];
            double rp = r_jph[j];

            double *piph = theDomain->piph[jk];
            double *dphi = theDomain->dphi[jk];
            double *gradr = theDomain->gradr[jk];

            for(i=0; i<Np[jk]; i++)
            {
                double phip = piph[i];
                double phim = phip - dphi[i];

                double xp[3] = {rp, phip, zp};
                double xm[3] = {rp, phim, zm};
                double dAp = get_dA(xp, xm, 1);
                xp[0] = rm;
                xm[0] = rm;
                double dAm = get_dA(xp, xm, 1);
                double idAtot = 1.0/(dAp + dAm);

                for(q=0; q<NUM_Q; q++)
                    gradr[NUM_Q*i+q] *= idAtot;
            }
        }
    }
    
    // Slope Limit
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr-1; j++)
        {
            int jkL = j+Nr*k;
            int jkR = j+1+Nr*k;
            int jkf = j + (Nr-1)*k;

            double *piphL = theDomain->piph[jkL];
            double *dphiL = theDomain->dphi[jkL];
            double *primL = theDomain->prim[jkL];
            double *gradrL = theDomain->gradr[jkL];
            double *gradpL = theDomain->gradp[jkL];
            
            double *piphR = theDomain->piph[jkR];
            double *dphiR = theDomain->dphi[jkR];
            double *primR = theDomain->prim[jkR];
            double *gradrR = theDomain->gradr[jkR];
            double *gradpR = theDomain->gradp[jkR];
            
            double *fr_phif = theDomain->fr_phif[jkf];
            double *fr_phib = theDomain->fr_phib[jkf];

            double rL = get_centroid(r_jph[j], r_jph[j-1], 1);
            double rR = get_centroid(r_jph[j+1], r_jph[j], 1);
            double idr = 1.0/(rR - rL);

            int iL = theDomain->I0[jkL];
            int iR = theDomain->I0[jkR];

            for( i=0 ; i<Nfr[jkf] ; ++i )
            {
                double dpL = get_signed_dp(0.5*(fr_phif[i]-fr_phib[i]), 
                                           piphL[iL] - 0.5*dphiL[iL]);
                double dpR = get_signed_dp(0.5*(fr_phif[i]-fr_phib[i]), 
                                           piphR[iR] - 0.5*dphiR[iR]);

                for( q=0 ; q<NUM_Q ; ++q )
                {
                    double fL = primL[NUM_Q*iL+q] + dpL * gradpL[NUM_Q*iL+q];
                    double fR = primR[NUM_Q*iR+q] + dpR * gradpR[NUM_Q*iR+q];
                    double gr = (fR - fL) * idr;
                    double grL = gradrL[NUM_Q*iL+q];
                    double grR = gradrR[NUM_Q*iR+q];
                    if(gr*grL < 0.0)
                        gradrL[NUM_Q*iL+q] = 0.0;
                    else if(fabs(PLM*gr) < fabs(grL))
                        gradrL[NUM_Q*iL+q] = PLM*gr;
                    if(gr*grR < 0.0)
                        gradrR[NUM_Q*iR+q] = 0.0;
                    else if(fabs(PLM*gr) < fabs(grR))
                        gradrR[NUM_Q*iR+q] = PLM*gr;
                }

                double dpLR = get_signed_dp(piphL[iL], piphR[iR]);
                if(dpLR < 0)
                    iL = iL==Np[jkL]-1 ? 0 : iL+1;
                else
                    iR = iR==Np[jkR]-1 ? 0 : iR+1;
            }
        }
}
