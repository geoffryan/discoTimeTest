
#include "header.h"
#include "domain_aos.h"
#include "geometry.h"
#include "hydro.h"
#include "faces.h"


void calc_dp_aos(struct domain *theDomain)
{
    struct cell ** theCells = theDomain->theCells;
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;

    int i,jk;
    for( jk=0 ; jk<Nr*Nz ; ++jk )
    {
        for( i=0 ; i<Np[jk] ; ++i )
        {
            int im = i-1;
            if( i == 0 )
                im = Np[jk]-1;
            double phim = theCells[jk][im].piph;
            double phip = theCells[jk][i ].piph;
            double dphi = get_dp(phip,phim);
            theCells[jk][i].dphi = dphi;
        }
    }
}

void setupGrid_aos(struct domain * theDomain)
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
}


void set_wcell_aos(struct domain *theDomain)
{
    struct cell ** theCells = theDomain->theCells;
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;

    int i,j,k;
    for( k=0 ; k<Nz ; ++k )
    {
        for( j=0 ; j<Nr ; ++j )
        {
            int jk = j+Nr*k;
            double w = 0.0;
            for( i=0 ; i<Np[jk] ; ++i )
                w += theCells[jk][i].wiph; 
            w /= (double)Np[jk];
            for( i=0 ; i<Np[jk] ; ++i )
                theCells[jk][i].wiph = w; 
        }    
    } 
}

void setupDomain_aos( struct domain * theDomain )
{
    setupGrid_aos(theDomain);

    srand(314159);
    rand();

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;
    theDomain->theCells = (struct cell **) malloc(Nr*Nz*sizeof(struct cell *));
   
    int jk;
    for( jk=0 ; jk<Nr*Nz ; ++jk )
        theDomain->theCells[jk] = (struct cell *) malloc(Np[jk]
                                                         *sizeof(struct cell));

    //Setup independent of node layout: pick the right rand()'s
    double Pmax = theDomain->phi_max;
    int i, j, k;


    for( k=0 ; k<Nz ; ++k )
    {
        //DO the work
        for( j=0 ; j<Nr ; ++j )
        {
            jk = k*Nr + j;
            double p0 = Pmax*(double)rand()/(double)RAND_MAX;
            double dp = Pmax/(double)Np[jk];
            for( i=0 ; i<Np[jk] ; ++i )
            {
                double phi = p0+dp*(double)i;
                if( phi > Pmax ) phi -= Pmax;
                    theDomain->theCells[jk][i].piph = phi;
                    theDomain->theCells[jk][i].dphi = dp;
            }
        }
    }

    theDomain->current_step = 0;
    theDomain->final_step  = theDomain->theParList.nsteps;

    theDomain->theFaces_1 = NULL;
    theDomain->theFaces_2 = NULL;
    theDomain->N_ftracks_r = get_num_rzFaces( Nr , Nz , 1 );
    theDomain->N_ftracks_z = get_num_rzFaces( Nr , Nz , 2 );

    theDomain->fIndex_r = (int *) malloc((theDomain->N_ftracks_r+1)
                                         *sizeof(int) );
    theDomain->fIndex_z = (int *) malloc((theDomain->N_ftracks_z+1)
                                         *sizeof(int) );

    setupCells_aos(theDomain);
}

void setupCells_aos( struct domain * theDomain ){

    calc_dp_aos( theDomain );

    int i,j,k;
    struct cell ** theCells = theDomain->theCells;
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;
    double * r_jph = theDomain->r_jph;
    double * z_kph = theDomain->z_kph;

    //Null setup for all cells
    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Nr; j++)
        {
            int jk = j+Nr*k;
            for(i=0; i<Np[jk]; i++)
            {
                struct cell * c = &(theCells[jk][i]);
                c->wiph = 0.0; 
                c->real = 0;
            }
        } 
    }

    //Setup all cells.
    int q;
    for( k=0 ; k<Nz ; ++k )
    {
        double z = get_centroid( z_kph[k], z_kph[k-1], 2);
        for( j=0 ; j<Nr ; ++j )
        {
            double r = get_centroid( r_jph[j], r_jph[j-1], 1);
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i )
            {
                struct cell * c = &(theCells[jk][i]);
                double phip = c->piph;
                double phim = phip-c->dphi;
                c->wiph = 0.0; 
                double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
                double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
                double phi = c->piph-.5*c->dphi;
                double x[3] = {r, phi, z};
    
                double dV = get_dV(xp, xm);
                initial( c->prim , x );
                prim2cons( c->prim , c->cons , x , dV, xp, xm);
                for(q=0; q<NUM_Q; q++)
                {
                    c->gradr[q] = 0.0;
                    c->gradp[q] = 0.0;
                    c->gradz[q] = 0.0;
                }
            }
        }    
    }    

    for( k=0 ; k<Nz ; ++k)
    {
      double z = get_centroid( z_kph[k], z_kph[k-1], 2);
      for( j=0 ; j<Nr ; ++j ){
         double r = get_centroid( r_jph[j], r_jph[j-1], 1);
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV( xp , xm );
            double phi = c->piph-.5*c->dphi;
            double x[3] = {r, phi, z};
            cons2prim( c->cons , c->prim , x , dV, xp, xm);
         }
      }
   }

   set_wcell_aos( theDomain );
}

void freeDomain_aos( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      free( theDomain->theCells[jk] );
   }
   free( theDomain->theCells );
   free( theDomain->Np );
   theDomain->r_jph--;
   free( theDomain->r_jph );
   theDomain->z_kph--;
   free( theDomain->z_kph );
   free( theDomain->fIndex_r );
   free( theDomain->fIndex_z );

   if(theDomain->theFaces_1 != NULL)
       free(theDomain->theFaces_1);
   if(theDomain->theFaces_2 != NULL)
       free(theDomain->theFaces_2);

}

double hash_aos(struct domain *theDomain, int qqq)
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
                    sum += theDomain->theCells[jk][i].prim[qqq]; 
                else
                {
                    int q;
                    for(q=0; q<NUM_Q; q++)
                        sum += theDomain->theCells[jk][i].prim[q]; 
                }
            }
        }

    return sum;
}

