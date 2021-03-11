#include "header.h"
#include "profiler.h"
#include "substep_aos.h"
#include "substep_soa1.h"
#include "timestep.h"

void timestep(struct domain *theDomain)
{
    prof_tick(theDomain->prof, PROF_RECON); 
#if TYPE == 0
    recon_aos(theDomain);
#elif TYPE == 1
    recon_soa1(theDomain);
#endif
    prof_tock(theDomain->prof, PROF_RECON); 


    prof_tick(theDomain->prof, PROF_FLUX); 
#if TYPE == 0
    addFlux_aos(theDomain);
#elif TYPE == 1
    addFlux_soa1(theDomain);
#endif
    prof_tock(theDomain->prof, PROF_FLUX); 
   

    prof_tick(theDomain->prof, PROF_SOURCE); 
#if TYPE == 0
    addSource_aos(theDomain);
#elif TYPE == 1
    addSource_soa1(theDomain);
#endif
    prof_tock(theDomain->prof, PROF_SOURCE); 

    prof_tick(theDomain->prof, PROF_C2P); 
#if TYPE == 0
    calcPrim_aos(theDomain);
#elif TYPE == 1
    calcPrim_soa1(theDomain);
#endif
    prof_tock(theDomain->prof, PROF_C2P); 

    theDomain->current_step += 1;
}
