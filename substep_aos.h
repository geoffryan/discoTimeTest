#ifndef DISCO_SUBSTEP_AOS_H
#define DISCO_SUBSTEP_AOS_H

#include "header.h"

void recon_aos(struct domain *theDomain);
void addFlux_aos(struct domain *theDomain, double dt);
void addSource_aos(struct domain *theDomain, double dt);
void calcPrim_aos(struct domain *theDomain);
void flux_phi_aos(struct domain *theDomain, double dt);
void flux_r_aos(struct domain *theDomain, double dt);
void flux_z_aos(struct domain *theDomain, double dt);
void riemann_phi_aos(struct cell *cL, struct cell *cR, double dt,
                     double rm, double rp, double r, double zm, double zp,
                     double z);
void riemann_r_aos(struct face *f, double dt, double r, double zp, double zm);
void riemann_z_aos(struct face *f, double dt, double rp, double rm, double z);

#endif
