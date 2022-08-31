#ifndef DISCO_DOMAIN_SOA2_H
#define DISCO_DOMAIN_SOA2_H

#include "header.h"

void setupDomain_soa2(struct domain *theDomain);
void build_faces_soa2(struct domain *theDomain);
void freeDomain_soa2(struct domain *theDomain);
double hash_soa2(struct domain *theDomain, int qqq);
void dump_grid_soa2(struct domain *theDomain);

#endif
