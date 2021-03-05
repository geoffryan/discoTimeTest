
#include "header.h"
#include "geometry.h"


double minmod_aos(double a, double b, double c)
{
    double m = a;
    if( a*b < 0.0 )
        m = 0.0;
    if( fabs(b) < fabs(m) )
        m = b;
    if( b*c < 0.0 )
        m = 0.0;
    if( fabs(c) < fabs(m) )
        m = c;
   
    return(m);
}

void plm_phi_aos( struct domain * theDomain )
{
    struct cell ** theCells = theDomain->theCells;
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
            for( i=0 ; i<Np[jk] ; ++i )
            {
                int im = i-1;
                if( im == -1 )
                    im = Np[jk]-1;
                int ip = (i+1)%Np[jk];
                struct cell * c  = &(theCells[jk][i]);
                struct cell * cL = &(theCells[jk][im]);
                struct cell * cR = &(theCells[jk][ip]);
                double dpL = cL->dphi;
                double dpC = c->dphi;
                double dpR = cR->dphi;
                for( q=0 ; q<NUM_Q ; ++q )
                {
                    double pL = cL->prim[q];
                    double pC = c->prim[q];
                    double pR = cR->prim[q];
                    double sL = pC - pL;
                    sL /= .5*( dpC + dpL );
                    double sR = pR - pC;
                    sR /= .5*( dpR + dpC );
                    double sC = pR - pL;
                    sC /= .5*( dpL + dpR ) + dpC;
                    c->gradp[q] = minmod_aos( PLM*sL , sC , PLM*sR );
                }
            }
        }
    }
}


void plm_trans_aos(struct domain *theDomain, struct face *theFaces, int Nf,
                   int dim )
{
    struct cell ** theCells = theDomain->theCells;
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;
    double * r_jph = theDomain->r_jph;
    double * z_kph = theDomain->z_kph;
    double PLM = 1.5;
   
    int i,j,k,q;

    //Clear gradients
    for( j=0 ; j<Nr ; ++j )
    {
        for( k=0 ; k<Nz ; ++k )
        {
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i )
            {
                if(dim == 1)
                    memset( theCells[jk][i].gradr , 0 , NUM_Q*sizeof(double) );
                else
                    memset( theCells[jk][i].gradz , 0 , NUM_Q*sizeof(double) );
            }
        }
    }

    if(PLM == 0.0)
        return;

    //Add weighted slopes
    int n;
    for( n=0 ; n<Nf ; ++n )
    {
        struct face * f  = &( theFaces[n] );
        double phi = f->cm[1];
        struct cell * cL = f->L;
        struct cell * cR = f->R;
        double dxL = f->dxL;
        double dxR = f->dxR;
        double phiL = cL->piph - .5*cL->dphi;
        double phiR = cR->piph - .5*cR->dphi;
        double dpL = get_signed_dp(phi,phiL);
        double dpR = get_signed_dp(phiR,phi);
        double dA  = f->dA;

        double *gradL = dim==1 ? cL->gradr : cL->gradz;
        double *gradR = dim==1 ? cR->gradr : cR->gradz;

        for( q=0 ; q<NUM_Q ; ++q )
        {
            double WL = cL->prim[q] + dpL*cL->gradp[q];
            double WR = cR->prim[q] - dpR*cR->gradp[q];

            double S = (WR-WL)/(dxR+dxL);
            gradL[q] += S*dA;
            gradR[q] += S*dA;
        }
    }

    //Divide by total weight
    for( j=0 ; j<Nr ; ++j )
    {
        for( k=0 ; k<Nz ; ++k )
        {
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i )
            {
                struct cell * c = &(theCells[jk][i]);
                double phip = c->piph;
                double phim = phip - c->dphi;
                double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
                double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
                if( dim==1 )
                    xm[0] = r_jph[j];
                else
                    xm[2] = z_kph[k];
                
                double dAp = get_dA(xp,xm,dim);
                if( dim==1 )
                {
                    xp[0] = r_jph[j-1];
                    xm[0] = r_jph[j-1];
                }
                else
                {
                    xp[2] = z_kph[k-1];
                    xm[2] = z_kph[k-1];
                }
                
                double dAm = get_dA(xp,xm,dim);
                double dAtot = dAp+dAm;
                if( (dim==1 && j==0   ) || (dim==2 && k==0   ) )
                    dAtot = dAp;
                if( (dim==1 && j==Nr-1) || (dim==2 && k==Nz-1) )
                    dAtot = dAm;
          
                double *grad = dim==1 ? c->gradr : c->gradz;
                for( q=0 ; q<NUM_Q ; ++q )
                   grad[q] /= dAtot;
            }    
        }    
    }

    //Slope Limiting
    for( n=0 ; n<Nf ; ++n )
    {
        struct face * f  = &( theFaces[n] );
        double phi = f->cm[1];
        struct cell * cL = f->L;
        struct cell * cR = f->R;
        double dxL = f->dxL;
        double dxR = f->dxR;
        double phiL = cL->piph - .5*cL->dphi;
        double phiR = cR->piph - .5*cR->dphi;
        double dpL = get_signed_dp(phi,phiL);
        double dpR = get_signed_dp(phiR,phi);
      
        double *gradL = dim==1 ? cL->gradr : cL->gradz;
        double *gradR = dim==1 ? cR->gradr : cR->gradz;

        for( q=0 ; q<NUM_Q ; ++q )
        {
            double WL = cL->prim[q] + dpL*cL->gradp[q];
            double WR = cR->prim[q] - dpR*cR->gradp[q];

            double S = (WR-WL)/(dxR+dxL);
            double SL = gradL[q];
            double SR = gradR[q];
            if( S*SL < 0.0 )
                gradL[q] = 0.0; 
            else if( fabs(PLM*S) < fabs(SL) )
                gradL[q] = PLM*S;
         
            if( S*SR < 0.0 )
                gradR[q] = 0.0; 
            else if( fabs(PLM*S) < fabs(SR) )
                gradR[q] = PLM*S;
         
        }    
    }
}

