#ifndef DISCO_SUBSTEP_SOA1_H
#define DISCO_SUBSTEP_SOA1_H

#include "header.h"

void recon_soa1(struct domain *theDomain);
void addFlux_soa1(struct domain *theDomain);
void addSource_soa1(struct domain *theDomain);
void calcPrim_soa1(struct domain *theDomain);

void flux_phi_soa1(struct domain *theDomain);
void flux_r_soa1(struct domain *theDomain);

void riemann_phi_soa1(struct domain *theDomain, int jk, int iL, int iR,
                      double dt, double rm, double rp, double r,
                      double zm, double zp, double z);
void riemann_phi_soa1_alt(double *prim, double *cons, double *dphi,
                      double *piph,
                      double *gradr, double *gradp, double *gradz,
                      int iL, int iR,
                      double dt, double rm, double rp, double r,
                      double zm, double zp, double z);
void riemann_r_soa1(struct domain *theDomain, int jkL, int jkR, int jkf,
                    int iL, int iR, int i, double dt, double rL, double rR,
                    double r, double zm, double zp, double z);
void riemann_r_soa1_alt(double *primL, double *primR,
                        double *consL, double *consR,
                        double *gradrL, double *gradpL, double *gradzL,
                        double *gradrR, double *gradpR, double *gradzR,
                        double *dphiL, double *piphL,
                        double *dphiR, double *piphR,
                        double *fr_dA, double *fr_phif, double *fr_phib,
                        int iL, int iR, int i, double dt, double rL, double rR,
                        double r, double zm, double zp, double z);
#endif
