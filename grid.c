#include "header.h"
#include "geometry.h"
#include "grid.h"

void setupGrid(struct domain * theDomain)
{

    int Ng = NUM_G;
    theDomain->Ng = Ng;
    int Num_R = theDomain->theParList.Num_R;
    int Num_Z = theDomain->theParList.Num_Z;

    double Rmin = theDomain->theParList.rmin;
    double Rmax = theDomain->theParList.rmax;
    double Zmin = theDomain->theParList.zmin;
    double Zmax = theDomain->theParList.zmax;
    double Pmax = theDomain->theParList.phimax;

    int N0r = 0;
    int N1r = Num_R;
    int NgRa = Ng;
    int NgRb = Ng;
    N0r -= NgRa;
    N1r += NgRb;
    int Nr = N1r-N0r;
    if(Num_R == 1)
    {
        NgRa = 0;
        NgRb = 0;
        N0r = 0;
        N1r = 1;
        Nr = 1;
    }

    int N0z = 0;
    int N1z = Num_Z;
    int NgZa = Ng;
    int NgZb = Ng;
    N0z -= NgZa;
    N1z += NgZb;
    int Nz = N1z-N0z;
    if(Num_Z == 1)
    {
        NgZa = 0;
        NgZb = 0;
        N0z = 0;
        N1z = 1;
        Nz = 1;
    }

    theDomain->Nr = Nr;
    theDomain->Nz = Nz;
    theDomain->NgRa = NgRa;
    theDomain->NgRb = NgRb;
    theDomain->NgZa = NgZa;
    theDomain->NgZb = NgZb;
    theDomain->N0r = N0r;
    theDomain->N0z = N0z;
    printf("Nr = %d, Nz = %d\n", Nr, Nz);

    theDomain->Np    = (int *)    malloc( Nr*Nz*sizeof(int) );
    theDomain->r_jph = (double *) malloc( (Nr+1)*sizeof(double) );
    theDomain->z_kph = (double *) malloc( (Nz+1)*sizeof(double) );

    ++(theDomain->r_jph);
    ++(theDomain->z_kph);

    int j,k;

    for( j=-1 ; j<Nr ; ++j )
    {
        double x = (N0r + j + 1) / (double) Num_R;
        theDomain->r_jph[j] = Rmin + x*(Rmax-Rmin);
    }
    
    double dz = (Zmax-Zmin)/(double)Num_Z;
    double z0 = Zmin + (double)N0z*dz;
    for( k=-1 ; k<Nz ; ++k )
        theDomain->z_kph[k] = z0 + ((double)k+1.)*dz;

    theDomain->phi_max = theDomain->theParList.phimax;
    setGeometryParams( theDomain );

    for(k=0 ; k<Nz ; ++k ){
        double zp = theDomain->z_kph[k];
        double zm = theDomain->z_kph[k-1];

        for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            double rp = theDomain->r_jph[j];
            double rm = theDomain->r_jph[j-1];

            double xp[3] = {rp, 0.0, zp};
            double xm[3] = {rm, 0.0, zm};
            double x[3];
            get_centroid_arr(xp, xm, x);

            double dr = get_dL(xp, xm, 1);
            double hp = get_scale_factor(xp, 0);

            double dp = dr/hp;
         
            int Np = (int)(Pmax/dp);
         
            if( Np<4 )
                Np=4;
            theDomain->Np[jk] = Np;
      }
   }

    theDomain->current_step = 0;
    theDomain->final_step  = theDomain->theParList.nsteps;
}

