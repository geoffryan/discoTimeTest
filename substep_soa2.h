#ifndef DISCO_SUBSTEP_SOA2_H
#define DISCO_SUBSTEP_SOA2_H

#include "header.h"

void recon_soa2(struct domain *theDomain);
void addFlux_soa2(struct domain *theDomain, double dt);
void addSource_soa2(struct domain *theDomain, double dt);
void calcPrim_soa2(struct domain *theDomain);

void flux_phi_soa2(struct domain *theDomain, double dt);
void flux_r_soa2(struct domain *theDomain, double dt);

void riemann_phi_soa2(struct domain *theDomain, int jk, int iL, int iR,
                      double dt, double rm, double rp, double r,
                      double zm, double zp, double z);
void riemann_soa2(const double *primL, const double *primR,
                  double *consL, double *consR,
                  const double *gradrL, const double *gradpL,
                  const double *gradzL, const double *gradrR,
                  const double *gradpR, const double *gradzR,
                  const double *xL, const double *xR,
                  double dt, double dA, const double *xp, const double *xm,
                  int dim);
#endif
