#ifndef DISCO_FACES_H
#define DISCO_FACES_H

#include "header.h"

int get_num_rzFaces(int Nr, int Nz, int dim);

void addFace(struct face *theFaces, int n, struct cell *cL, struct cell *cR,
             double dxL, double dxR, double *xp, double *xm, int dim,
             int LRtype);
void buildfaces(struct domain *theDomain, int dim, int mode);
void setup_faces(struct domain *theDomain, int dim);

#endif
