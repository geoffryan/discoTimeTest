#ifndef DISCO_PLMAOS_H
#define DISCO_PLMAOS_H

#include "header.h"

double minmod_aos(double a, double b, double c);

void plm_phi_aos( struct domain * theDomain );


void plm_trans_aos(struct domain *theDomain, struct face *theFaces, int Nf,
                   int dim );

#endif
