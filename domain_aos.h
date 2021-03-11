#ifndef DISCO_DOMAIN_H
#define DISCO_DOMAIN_H

#include "header.h"

void calc_dp_aos(struct domain *);
void set_wcell_aos(struct domain *);
void setupGrid_aos(struct domain *theDomain);
void setupDomain_aos(struct domain *theDomain);
void setupCells_aos(struct domain *theDomain);
void freeDomain_aos(struct domain *theDomain);
double hash_aos(struct domain *theDomain, int qqq);
void dump_grid_aos(struct domain *theDomain);

#endif
