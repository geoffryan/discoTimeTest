#ifndef DISCO_DOMAIN_SOA2_H
#define DISCO_DOMAIN_SOA2_H

#include "header.h"

void setupDomain_soa2(struct domain *theDomain);
void build_faces_soa2(struct domain *theDomain, int dim);
void freeDomain_soa2(struct domain *theDomain);
double hash_soa2(struct domain *theDomain, int qqq);
void dump_grid_soa2(struct domain *theDomain, char label[]);

void fill_faces_soa2(int Np, int Np1, const double *piph, const double *piph1,
                     const double *dphi, const double *dphi1,
                     double rm, double rp, double zm, double zp,
                     int *f_i, double *f_phib, double *f_dA, double *f_dphi,
                     double phi_max, int inc_match, int dim);
void check_faces(struct domain *theDomain);

#endif
