#include "header.h"
#include "profiler.h"
#include "domain_aos.h"
#include "domain_soa1.h"
#include "domain_soa2.h"
#include "substep_aos.h"
#include "substep_soa1.h"
#include "substep_soa2.h"
#include "timestep.h"

void dump(struct domain *theDomain, char label[])
{
#if DEBUG
#if TYPE == 0
    dump_grid_aos(theDomain, label);
#elif TYPE == 1
    dump_grid_soa1(theDomain, label);
#elif TYPE == 2
    dump_grid_soa2(theDomain, label);
#endif
#endif
}

void timestep(struct domain *theDomain)
{
    prof_tick(theDomain->prof, PROF_RECON); 
#if TYPE == 0
    recon_aos(theDomain);
#elif TYPE == 1
    recon_soa1(theDomain);
#elif TYPE == 2
    recon_soa2(theDomain);
#endif
    prof_tock(theDomain->prof, PROF_RECON);
    dump(theDomain, "recons");

    double dt = theDomain->theParList.dt;


    prof_tick(theDomain->prof, PROF_FLUX); 
#if TYPE == 0
    addFlux_aos(theDomain, dt);
#elif TYPE == 1
    addFlux_soa1(theDomain, dt);
#elif TYPE == 2
    addFlux_soa2(theDomain, dt);
#endif
    prof_tock(theDomain->prof, PROF_FLUX); 
    dump(theDomain, "fluxes");
   

    prof_tick(theDomain->prof, PROF_SOURCE); 
#if TYPE == 0
    addSource_aos(theDomain, dt);
#elif TYPE == 1
    addSource_soa1(theDomain, dt);
#elif TYPE == 2
    addSource_soa2(theDomain, dt);
#endif
    prof_tock(theDomain->prof, PROF_SOURCE); 
    dump(theDomain, "source");

    prof_tick(theDomain->prof, PROF_C2P); 
#if TYPE == 0
    calcPrim_aos(theDomain);
#elif TYPE == 1
    calcPrim_soa1(theDomain);
#elif TYPE == 2
    calcPrim_soa2(theDomain);
#endif
    prof_tock(theDomain->prof, PROF_C2P); 
    dump(theDomain, "con2pr");

    theDomain->current_step += 1;
}
