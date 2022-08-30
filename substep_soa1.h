#ifndef DISCO_SUBSTEP_SOA1_H
#define DISCO_SUBSTEP_SOA1_H

#include "header.h"

void recon_soa1(struct domain *theDomain);
void addFlux_soa1(struct domain *theDomain, double dt);
void addSource_soa1(struct domain *theDomain, double dt);
void calcPrim_soa1(struct domain *theDomain);

void flux_phi_soa1(struct domain *theDomain, double dt);
void flux_r_soa1(struct domain *theDomain, double dt);

void riemann_phi_soa1(struct domain *theDomain, int jk, int iL, int iR,
                      double dt, double rm, double rp, double r,
                      double zm, double zp, double z);
void riemann_phi_soa1_alt(const double *prim, double *cons,
                          const double *dphi, const double *piph,
                          const double *gradr, const double *gradp,
                          const double *gradz,
                          int iL, int iR,
                          double dt, double rm, double rp, double r,
                          double zm, double zp, double z);
void riemann_r_soa1(struct domain *theDomain, int jkL, int jkR, int jkf,
                    int iL, int iR, int i, double dt, double rL, double rR,
                    double r, double zm, double zp, double z);
void riemann_r_soa1_alt(const double *primL, const double *primR,
                        double *consL, double *consR,
                        const double *gradrL, const double *gradpL,
                        const double *gradzL,
                        const double *gradrR, const double *gradpR,
                        const double *gradzR,
                        const double *dphiL, const double *piphL,
                        const double *dphiR, const double *piphR,
                        const double *fr_dA, const double *fr_phif,
                        const double *fr_phib,
                        int iL, int iR, int i, double dt, double rL, double rR,
                        double r, double zm, double zp, double z);
#endif
