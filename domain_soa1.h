#ifndef DISCO_DOMAIN_SOA1_H
#define DISCO_DOMAIN_SOA1_H

#include "header.h"

void setupDomain_soa1(struct domain *theDomain);
void build_faces_soa1(struct domain *theDomain);
void freeDomain_soa1(struct domain *theDomain);
double hash_soa1(struct domain *theDomain, int qqq);
void dump_grid_soa1(struct domain *theDomain, char label[]);

#endif
