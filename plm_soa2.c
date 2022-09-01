
#include "header.h"
#include "geometry.h"
#include "plm_soa2.h"

double minmod_soa2(double a, double b, double c)
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

void plm_phi_soa2(struct domain *theDomain)
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
                    gradp[iqC+q] = minmod_soa2( PLM*sL , sC , PLM*sR );
                }
            }
        }
    }
}

void plm_r_soa2(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

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
        for(j=0; j<Nr; j++)
        {
            int jk = j+Nr*k;
            int jkL = jk-1;
            int jkR = jk+1;

            double *piph = theDomain->piph[jk];
            double *dphi = theDomain->dphi[jk];
            double *prim = theDomain->prim[jk];
            double *gradr = theDomain->gradr[jk];
            double *gradp = theDomain->gradp[jk];

            double r = get_centroid(r_jph[j], r_jph[j-1], DIM_R);

            if(j > 0)
            {
                double rL = get_centroid(r_jph[j-1], r_jph[j-2], DIM_R);

                double *dphiL = theDomain->dphi[jkL];
                double *primL = theDomain->prim[jkL];
                double *gradrL = theDomain->gradr[jkL];
                double *gradpL = theDomain->gradp[jkL];

                int *fr_iL = theDomain->fr_iL[jk];
                double *fr_dA_L = theDomain->fr_dA_L[jk];
                double *fr_dphi_L = theDomain->fr_dphi_L[jk];
                double *fr_offset_L = theDomain->fr_offset_L[jk];
                    
                double idr = 1.0 / (r - rL);

                for(i=0; i<Np[jk]; i++)
                {
                    int iL = fr_iL[i];
                    double phi = piph[i] - 0.5*dphi[i];
                    double phif = piph[i] - 0.5*fr_dphi_L[i];
                    double phiL = piph[i] + fr_offset_L[i] - 0.5*dphiL[iL];

                    double dpL = phif - phiL;
                    double dpR = phif - phi;

                    for(q=0; q<NUM_Q; q++)
                    {
                        int iq = NUM_Q*i + q;
                        int iqL = NUM_Q*iL + q;
                        double fL = primL[iqL] + dpL * gradpL[iqL];
                        double fR = prim[iq] + dpR * gradp[iq];
                        double gr = (fR - fL) * idr;
                        gradrL[iqL] += gr * fr_dA_L[i];
                        gradr[iq] += gr * fr_dA_L[i];
                    }
                }
            }

            if(j < Nr-1)
            {
                double rR = get_centroid(r_jph[j+1], r_jph[j], DIM_R);

                double *dphiR = theDomain->dphi[jkR];
                double *primR = theDomain->prim[jkR];
                double *gradrR = theDomain->gradr[jkR];
                double *gradpR = theDomain->gradp[jkR];

                int *fr_iR = theDomain->fr_iR[jk];
                double *fr_dA_R = theDomain->fr_dA_R[jk];
                double *fr_dphi_R = theDomain->fr_dphi_R[jk];
                double *fr_offset_R = theDomain->fr_offset_R[jk];
                    
                double idr = 1.0 / (rR - r);

                for(i=0; i<Np[jk]; i++)
                {
                    int iR = fr_iR[i];
                    double phi = piph[i] - 0.5*dphi[i];
                    double phif = piph[i] - 0.5*fr_dphi_R[i];
                    double phiR = piph[i] + fr_offset_R[i] - 0.5*dphiR[iR];

                    double dpL = phif - phi;
                    double dpR = phif - phiR;

                    for(q=0; q<NUM_Q; q++)
                    {
                        int iq = NUM_Q*i + q;
                        int iqR = NUM_Q*iR + q;
                        double fL = prim[iq] + dpL * gradp[iq];
                        double fR = primR[iqR] + dpR * gradpR[iqR];
                        double gr = (fR - fL) * idr;
                        gradr[iq] += gr * fr_dA_R[i];
                        gradrR[iqR] += gr * fr_dA_R[i];
                    }
                }
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
                double idAtot;
                if(j == 0)
                    idAtot = 1.0/dAp;
                else if(j == Nr-1)
                    idAtot = 1.0/dAm;
                else
                    idAtot = 1.0/(dAp + dAm);

                for(q=0; q<NUM_Q; q++)
                    gradr[NUM_Q*i+q] *= idAtot;
            }
        }
    }
    
    // Slope Limit
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = j+Nr*k;
            int jkL = jk-1;
            int jkR = jk+1;

            double *piph = theDomain->piph[jk];
            double *dphi = theDomain->dphi[jk];
            double *prim = theDomain->prim[jk];
            double *gradr = theDomain->gradr[jk];
            double *gradp = theDomain->gradp[jk];

            double r = get_centroid(r_jph[j], r_jph[j-1], DIM_R);

            if(j > 0)
            {
                double rL = get_centroid(r_jph[j-1], r_jph[j-2], DIM_R);

                double *dphiL = theDomain->dphi[jkL];
                double *primL = theDomain->prim[jkL];
                double *gradrL = theDomain->gradr[jkL];
                double *gradpL = theDomain->gradp[jkL];

                int *fr_iL = theDomain->fr_iL[jk];
                double *fr_dphi_L = theDomain->fr_dphi_L[jk];
                double *fr_offset_L = theDomain->fr_offset_L[jk];
                    
                double idr = 1.0 / (r - rL);

                for(i=0; i<Np[jk]; i++)
                {
                    int iL = fr_iL[i];
                    double phi = piph[i] - 0.5*dphi[i];
                    double phif = piph[i] - 0.5*fr_dphi_L[i];
                    double phiL = piph[i] + fr_offset_L[i] - 0.5*dphiL[iL];

                    double dpL = phif - phiL;
                    double dpR = phif - phi;

                    for(q=0; q<NUM_Q; q++)
                    {
                        int iq = NUM_Q*i + q;
                        int iqL = NUM_Q*iL + q;
                        double fL = primL[iqL] + dpL * gradpL[iqL];
                        double fR = prim[iq] + dpR * gradp[iq];
                        double gr = (fR - fL) * idr;
                        double grL = gradrL[iqL];
                        double grR = gradr[iq];
                        if(gr*grL < 0.0)
                            gradrL[iqL] = 0.0;
                        else if(fabs(PLM*gr) < fabs(grL))
                            gradrL[iqL] = PLM*gr;

                        if(gr*grR < 0.0)
                            gradr[iq] = 0.0;
                        else if(fabs(PLM*gr) < fabs(grR))
                            gradr[iq] = PLM*gr;
                    }
                }
            }

            if(j < Nr-1)
            {
                double rR = get_centroid(r_jph[j+1], r_jph[j], DIM_R);

                double *dphiR = theDomain->dphi[jkR];
                double *primR = theDomain->prim[jkR];
                double *gradrR = theDomain->gradr[jkR];
                double *gradpR = theDomain->gradp[jkR];

                int *fr_iR = theDomain->fr_iR[jk];
                double *fr_dphi_R = theDomain->fr_dphi_R[jk];
                double *fr_offset_R = theDomain->fr_offset_R[jk];
                    
                double idr = 1.0 / (rR - r);

                for(i=0; i<Np[jk]; i++)
                {
                    int iR = fr_iR[i];
                    double phi = piph[i] - 0.5*dphi[i];
                    double phif = piph[i] - 0.5*fr_dphi_R[i];
                    double phiR = piph[i] + fr_offset_R[i] - 0.5*dphiR[iR];

                    double dpL = phif - phi;
                    double dpR = phif - phiR;

                    for(q=0; q<NUM_Q; q++)
                    {
                        int iq = NUM_Q*i + q;
                        int iqR = NUM_Q*iR + q;
                        double fL = prim[iq] + dpL * gradp[iq];
                        double fR = primR[iqR] + dpR * gradpR[iqR];
                        double gr = (fR - fL) * idr;
                        double grL = gradr[iq];
                        double grR = gradrR[iqR];
                        if(gr*grL < 0.0)
                            gradr[iq] = 0.0;
                        else if(fabs(PLM*gr) < fabs(grL))
                            gradr[iq] = PLM*gr;

                        if(gr*grR < 0.0)
                            gradrR[iqR] = 0.0;
                        else if(fabs(PLM*gr) < fabs(grR))
                            gradrR[iqR] = PLM*gr;
                    }
                }
            }
        }
}
