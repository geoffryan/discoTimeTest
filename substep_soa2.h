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
void riemann_r_soa2(struct domain *theDomain, int jkL, int jkR, int jkf,
                    int iL, int iR, int i, double dt, double rL, double rR,
                    double r, double zm, double zp, double z);
#endif
