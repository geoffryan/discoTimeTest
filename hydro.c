
#include "header.h"
#include "hydro.h"
#include "geometry.h"

static double gamma_law = 5.0/3.0; 

int set_B_flag()
{
    return 0;
}

double get_omega(const double *prim, const double *x)
{
    return prim[UPP];
}

void initial(double *prim, const double *x)
{
    int q;
    double r = x[0];
    double phi = x[1];
    double z = x[3];
    prim[RHO] = 1.0;
    prim[PPP] = 1.0 + 0.01 * r * sin(4*phi);
    prim[URR] = r;
    prim[UZZ] = z;
    prim[UPP] = 0.0;
    if(NUM_N > 0)
        prim[NUM_C] = log(prim[PPP] * pow(prim[RHO], -gamma_law))
                         / (gamma_law - 1);
    
    for(q=NUM_C+1; q<NUM_Q; q++)
        prim[q] = q + r*cos(phi) + z*z*sin(phi);
}

void prim2cons(const double * prim, double * cons, const double * x,
               double dV, const double * xp, const double *xm){

    double r = x[0];
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vr  = prim[URR];
    double vp  = prim[UPP]*r;
    double vz  = prim[UZZ];

    double v2  = vr*vr + vp*vp + vz*vz;

    double rhoe = Pp/(gamma_law - 1.);

    cons[DDD] = rho*dV;
    cons[TAU] = (.5*rho*v2 + rhoe )*dV;
    cons[SRR] = rho*vr*dV;
    cons[LLL] = r*rho*vp*dV;
    cons[SZZ] = rho*vz*dV;

    int q;
    for( q=NUM_C ; q<NUM_Q ; ++q )
        cons[q] = prim[q]*cons[DDD];
}

void cons2prim(const double *cons, double *prim, const double *x, double dV,
               const double *xp, const double *xm)
{
   
    double r = x[0];

    double rho = cons[DDD]/dV;
    double Sr  = cons[SRR]/dV;
    double Sp  = cons[LLL]/(dV * r);
    double Sz  = cons[SZZ]/dV;
    double E   = cons[TAU]/dV;
   
    double vr = Sr/rho;
    double vp = Sp/rho;
    double vz = Sz/rho;

    double KE = .5*( Sr*vr + rho*vp*vp + Sz*vz );
    double rhoe = E-KE;
    double Pp = (gamma_law - 1.)*rhoe;

    if(NUM_N > 0)
    {
        //EXPERIMENTAL ENTROPY FIX
        double s = cons[NUM_C] / cons[DDD];
        Pp = exp((gamma_law-1.0)*s) * pow(rho, gamma_law);
    }

    prim[RHO] = rho;
    prim[PPP] = Pp;
    prim[URR] = vr;
    prim[UPP] = vp/r;
    prim[UZZ] = vz;

    int q;
    for( q=NUM_C ; q<NUM_Q ; ++q )
        prim[q] = cons[q]/cons[DDD];
}

void flux(const double * prim, double * flux, const double * x,
          const double * n, const double *xp, const double *xm)
{
    double r = x[0];
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vr  = prim[URR];
    double vp  = prim[UPP]*r;
    double vz  = prim[UZZ];

    double vn = vr*n[0] + vp*n[1] + vz*n[2];

    double rhoe = Pp/(gamma_law - 1.);
    double v2 = vr*vr + vp*vp + vz*vz;

    int q;
    flux[DDD] = rho*vn;
    flux[SRR] = rho*vr*vn + Pp*n[0];
    flux[LLL] = r*(rho*vp*vn + Pp*n[1]);
    flux[SZZ] = rho*vz*vn + Pp*n[2];
    flux[TAU] = ( .5*rho*v2 + rhoe + Pp )*vn - Pp*vn;

    for( q=NUM_C ; q<NUM_Q ; ++q )
        flux[q] = prim[q]*flux[DDD];
}

void source(const double *prim, double *cons, const double *xp,
            const double *xm, double dVdt)
{
    double rp = xp[0];
    double rm = xm[0];
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double r = get_centroid(rp, rm, 1);
    double r_1  = .5*(rp+rm);
    double omega = prim[UPP];

    double centrifugal = rho*omega*omega*r;

    double press_bal   = Pp/r_1;

    // Geometric source term
    cons[SRR] += dVdt*( centrifugal + press_bal );
}

void visc_flux(const double * prim, const double * gradr, const double * gradp,
               const double * gradz, double * flux,
               const double * x, const double * n)
{
   double r = x[0];
   double nu = 1.0e-6;

   double rho = prim[RHO];
   double vr  = prim[URR];
   double om  = prim[UPP];
   double vz  = prim[UZZ];

   //Divergence of v divided by number of spatial dimensions (3)
   double divV_o_d = (gradr[URR] + gradp[UPP] + gradz[UZZ] + vr/r) / 3.0;

   // Covariant components of shear tensor.
   double srr = gradr[URR] - divV_o_d;
   double spp = r*(r*gradp[UPP] + vr - r*divV_o_d);
   double szz = gradz[UZZ] - divV_o_d;
   double srp = 0.5*(r*r*gradr[UPP] + gradp[URR]);
   double srz = 0.5*(gradr[UZZ] + gradz[URR]);
   double spz = 0.5*(gradp[UZZ] + r*r*gradz[UPP]);

   // Covariant components of shear normal to face, shear_{ij} * n^{j}.
   // Given n is in orthonormal basis, 1/r factor corrects to coordinate basis
   double nc[3] = {n[0], n[1]/r, n[2]};
   double srn = srr*nc[0] + srp*nc[1] + srz*nc[2];
   double spn = srp*nc[0] + spp*nc[1] + spz*nc[2];
   double szn = srz*nc[0] + spz*nc[1] + szz*nc[2];

   flux[SRR] = -2 * nu * rho * srn;
   flux[LLL] = -2 * nu * rho * spn;
   flux[SZZ] = -2 * nu * rho * szn;
   flux[TAU] = -2 * nu * rho * ( vr*srn + om*spn + vz*szn);
}

void visc_source(const double * prim, const double * gradr, const double *gradp,
                 const double * gradz, double * cons, const double *xp,
                 const double *xm, double dVdt)
{
   double x[3];
   get_centroid_arr(xp, xm, x);
   double r = x[0];
   double nu = 1.0e-6;

   double rho = prim[RHO];
   double vr  = prim[URR];

   //Divergence of v divided by number of spatial dimensions (3)
   double divV_o_d = (gradr[URR] + gradp[UPP] + gradz[UZZ] + vr/r) / 3.0;

   // Contravariant -phi-phi component of shear tensor.
   double spp = (r*gradp[UPP] + vr - r*divV_o_d) / (r*r*r);
   
   cons[SRR] += (-2 * r * rho * nu * spp) * dVdt;

   // Mixed ^r_phi and ^z_phi components of shear tensor
   double srp = 0.5*(r*r*gradr[UPP] + gradp[URR]);
   double spz = 0.5*(gradp[UZZ] + r*r*gradz[UPP]);

   if(NUM_N > 0)
   {
       //EXPERIMENTAL Viscous entropy source
       //
       // covariant _{ij} shear tensor components that we don't
       // already have
       double sRR = gradr[URR] - divV_o_d;
       double sZZ = gradz[UZZ] - divV_o_d;
       double sRZ = 0.5*(gradr[UZZ] + gradz[URR]);

       double s2 = sRR*sRR + spp*spp * (r*r*r*r) + sZZ*sZZ
                    + 2*(srp*srp + spz*spz)/(r*r) + 2*sRZ*sRZ;

       // IDEAL GAS TEMPERATURE. NOT APPROPRIATE FOR RAD-PRESSURE GAS
       double T = prim[PPP] / rho;

       cons[NUM_C] += 2*rho*nu * s2 / T * dVdt;
   }
}

void flux_to_E( const double * Flux , const double * Ustr , const double * x, 
                double * E1_riemann , double * B1_riemann , 
                double * E2_riemann , double * B2_riemann , int dim ){

   //Silence is Golden.

}

void vel( const double * prim1 , const double * prim2 , 
         double * Sl , double * Sr , double * Ss , 
         const double * n , const double * x , double * Bpack ){

    double r = x[0];
    double P1   = prim1[PPP];
    double rho1 = prim1[RHO];
    double vn1  = prim1[URR]*n[0] + prim1[UPP]*n[1]*r + prim1[UZZ]*n[2];

    double cs1 = sqrt(gamma_law*P1/rho1);

    double P2   = prim2[PPP];
    double rho2 = prim2[RHO];
    double vn2  = prim2[URR]*n[0] + prim2[UPP]*n[1]*r + prim2[UZZ]*n[2];

    double cs2 = sqrt(gamma_law*P2/rho2);

    *Ss = (P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )
           /( rho1*(-cs1) - rho2*cs2 );

    *Sr =  cs1 + vn1;
    *Sl = -cs1 + vn1;

    if(*Sr <  cs2 + vn2 )
        *Sr =  cs2 + vn2;
    if(*Sl > -cs2 + vn2 )
        *Sl = -cs2 + vn2;

}

double mindt(const double * prim , double w ,
             const double * xp , const double * xm ){

   double r = get_centroid(xp[0], xm[0], 1);
   double Pp  = prim[PPP];
   double rho = prim[RHO];
   double vp  = (prim[UPP]-w)*r;
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double cs  = sqrt(gamma_law*Pp/rho);

   double maxvr = cs + fabs(vr);
   double maxvp = cs + fabs(vp);
   double maxvz = cs + fabs(vz);

   double dtr = get_dL(xp,xm,1)/maxvr;
   double dtp = get_dL(xp,xm,0)/maxvp;
   double dtz = get_dL(xp,xm,2)/maxvz;
   
   double dt = dtr;
   if( dt > dtp ) dt = dtp;
   if( dt > dtz ) dt = dtz;

   return( dt );

}

void reflect_prims(double * prim, const double * x, int dim)
{
    //dim == 0: r, dim == 1: p, dim == 2: z
    if(dim == 0)
        prim[URR] = -prim[URR];
    else if(dim == 1)
        prim[UPP] = -prim[UPP];
    else if(dim == 2)
        prim[UZZ] = -prim[UZZ];
}

double bfield_scale_factor(double x, int dim)
{
    // Returns the factor used to scale B_cons.
    // x is coordinate location in direction dim.
    // dim == 0: r, dim == 1: p, dim == 2: z
    
    return 1.0;
}
