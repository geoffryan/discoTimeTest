#ifndef DISCO_GEOMETRY_H
#define DISCO_GEOMETRY_H

#include "header.h"

void setGeometryParams( struct domain * theDomain );

double get_dp( double phip , double phim );
double get_signed_dp( double phip , double phim );

double get_centroid(double xp, double xm, int dim);
double get_dL( const double * xp , const double * xm , int dim );
double get_dA( const double * xp , const double * xm , int dim );
double get_dV( const double * xp , const double * xm );
double get_scale_factor( const double * x, int dim);
double get_vol_element(const double *x);

void get_xyz(const double *x, double *xyz);
void get_rpz(const double *x, double *rpz);
void get_coords_from_xyz(const double *xyz, double *x);
void get_coords_from_rpz(const double *rpz, double *x);

void get_vec_rpz(const double *x, const double *v, double *vrpz);
void get_vec_from_rpz(const double *x, const double *vrpz, double *v);
void get_vec_xyz(const double *x, const double *v, double *vxyz);
void get_vec_from_xyz(const double *x, const double *vxyz, double *v);

void geom_grad(const double *prim, double *grad,
               const double *xp,const double *xm, double PLM, int dim, int LR);
void geom_interpolate(const double *prim, const double *gradp,
                      const double *gradT, const double *x,
                      double dphi, double dxT, double * primI, int dim);

void get_centroid_arr(const double *xp, const double *xm, double *x);
void get_vec_contravariant(const double *x, double *v, double *vc);
void get_vec_covariant(const double *x, double *v, double *vc);

void geom_polar_vec_adjust(const double *xp, const double *xm, double *fac);

#endif
