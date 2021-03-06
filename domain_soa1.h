#ifndef DISCO_DOMAIN_SOA1_H
#define DISCO_DOMAIN_SOA1_H

#include "header.h"

void setupDomain_soa1(struct domain *theDomain);
void freeDomain_soa1(struct domain *theDomain);
double hash_soa1(struct domain *theDomain, int qqq);

#endif
