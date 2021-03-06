
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
