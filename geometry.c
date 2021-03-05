
#include "header.h"

static double phi_max = 0.0;

void setGeometryParams( struct domain * theDomain ){
   phi_max = theDomain->phi_max;
}

double get_dp( double phip , double phim ){
   double dp = phip-phim;
   while( dp < 0.0 ) dp += phi_max;
   while( dp > phi_max) dp -= phi_max;
   return(dp);
}

double get_signed_dp( double phip , double phim ){
   double dp = phip-phim;
   while( dp <-.5*phi_max ) dp += phi_max;
   while( dp > .5*phi_max ) dp -= phi_max;
   return(dp);
}

double get_centroid(double xp, double xm, int dim)
{
    if(dim == 1)
        return 2.0*(xp*xp + xp*xm + xm*xm) / (3.0*(xp+xm));
    else
        return 0.5*(xp+xm);
}

double get_dL( const double * xp , const double * xm , int dim ){
    double r = .5*(xp[0]+xm[0]);
    double dphi = get_dp(xp[1], xm[1]);
    if(dim == 0)
        return r*dphi;
    else if(dim == 1)
        return xp[0]-xm[0];
    else
        return xp[2]-xm[2];
}

double get_dA( const double * xp , const double * xm , int dim ){
    double r  = .5*(xp[0]+xm[0]);
    double dr   = xp[0]-xm[0];
    double dphi = get_dp(xp[1], xm[1]);
    double dz   = xp[2]-xm[2];

    if(dim == 0)
        return dr*dz;
    else if(dim == 1)
        return r*dphi*dz;
    else
        return r*dr*dphi;
}

double get_dV( const double * xp , const double * xm ){
    double r  = .5*(xp[0]+xm[0]);
    double dr   = xp[0]-xm[0];
    double dphi = get_dp(xp[1],xm[1]);
    double dz   = xp[2]-xm[2];

    return( r*dr*dphi*dz );
}

double get_scale_factor( const double * x, int dim)
{
    if(dim == 0)
        return x[0];
    return 1.0;
}

double get_vol_element(const double *x)
{
    return x[0];
}

void get_xyz(const double *x, double *xyz)
{
    xyz[0] = x[0] * cos(x[1]);
    xyz[1] = x[0] * sin(x[1]);
    xyz[2] = x[2];
}

void get_rpz(const double *x, double *rpz)
{
    rpz[0] = x[0];
    rpz[1] = x[1];
    rpz[2] = x[2];
}

void get_coords_from_xyz(const double *xyz, double *x)
{
    x[0] = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
    x[1] = atan2(xyz[1],xyz[0]);
    x[2] = xyz[2];
}

void get_coords_from_rpz(const double *rpz, double *x)
{
    x[0] = rpz[0];
    x[1] = rpz[1];
    x[2] = rpz[2];
}

void get_vec_rpz(const double *x, const double *v, double *vrpz)
{
    vrpz[0] = v[0];
    vrpz[1] = v[1];
    vrpz[2] = v[2];
}

void get_vec_from_rpz(const double *x, const double *vrpz, double *v)
{
    v[0] = vrpz[0];
    v[1] = vrpz[1];
    v[2] = vrpz[2];
}

void get_vec_xyz(const double *x, const double *v, double *vxyz)
{
    double phi = x[1];
    double cp = cos(phi);
    double sp = sin(phi);

    vxyz[0] = cp*v[0] - sp*v[1];
    vxyz[1] = sp*v[0] + cp*v[1];
    vxyz[2] = v[2];
}

void get_vec_from_xyz(const double *x, const double *vxyz, double *v)
{
    double phi = x[1];
    double cp = cos(phi);
    double sp = sin(phi);

    v[0] =  cp*vxyz[0] + sp*vxyz[1];
    v[1] = -sp*vxyz[0] + cp*vxyz[1];
    v[2] = vxyz[2];
}

void geom_grad(const double *prim, double *grad, const double *xp, const double *xm, 
                double PLM, int dim, int LR)
{
    if(dim !=1 || LR != 0)
    {
        printf("Geometric gradient called on non-geometric boundary\n");
        printf("--Cylindrical setup only has geometric boundary at r=0.\n");
        return;
    }
    if(xp[0] < 0.0 || fabs(xm[0]) > 1.0e-10*xp[0])
    {
        printf("Geometric gradient called on cell with rm = %.le (!= 0)\n",
                xm[0]);
        return;
    }

    int q;
    double r = get_centroid(xp[0], xm[0], 1);
    for(q = 0; q<NUM_Q; q++)
    {
        if(q == URR || (NUM_C>BZZ && q==BRR))
        {
            double SL = prim[q]/r;
            double S = grad[q];
            if( S*SL < 0.0 )
                grad[q] = 0.0; 
            else if( fabs(PLM*SL) < fabs(S) )
                grad[q] = PLM*SL;
        }
        else
            grad[q] = 0.0;
    }
}

void geom_polar_vec_adjust(const double *xp, const double *xm, double *fac)
{
    double dphi = get_dp(xp[1], xm[1]);
    double adjust = sin(0.5*dphi) / (0.5*dphi);

    fac[0] = adjust;
    fac[1] = adjust;
    fac[2] = 1.0;
}

void geom_interpolate(const double *prim, const double *gradp,
                      const double *gradT, const double *x,
                      double dphi, double dxT, double * primI, int dim)
{
    double r = x[0];
    double sp = sin(dphi);
    double cp = cos(dphi);

    primI[RHO] = prim[RHO] + dphi * gradp[RHO];
    primI[PPP] = prim[PPP] + dphi * gradp[PPP];

    primI[URR] = cp*prim[URR] + sp*r*prim[UPP];
    primI[UPP] = -sp*prim[URR]/r + cp*prim[UPP];

    primI[UZZ] = prim[UZZ] + dphi * gradp[UZZ];

    int q;
    for(q=5; q<NUM_Q; q++)
        primI[q] = prim[q] + dphi * gradp[q];

    if(dim == 1)
    {
        primI[RHO] = prim[RHO] + dxT * gradT[RHO];
        primI[PPP] = prim[PPP] + dxT * gradT[PPP];

        primI[UPP] = primI[UPP]*r/(r+dxT);

        primI[UZZ] = prim[UZZ] + dxT * gradT[UZZ];

        for(q=5; q<NUM_Q; q++)
            primI[q] = prim[q] + dphi * gradp[q];
    }

    if(dim == 2)
    {
        primI[RHO] = prim[RHO] + dxT * gradT[RHO];
        primI[PPP] = prim[PPP] + dxT * gradT[PPP];

        primI[UZZ] = prim[UZZ] + dxT * gradT[UZZ];

        for(q=5; q<NUM_Q; q++)
            primI[q] = prim[q] + dphi * gradp[q];
    }
}
