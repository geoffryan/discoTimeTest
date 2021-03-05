
#include "header.h"
#include "domain_aos.h"
#include "profiler.h"
#include "timestep.h"

void read_par_file( struct domain * );

int main( int argc , char * argv[] )
{
 
    struct domain theDomain = {0};
    struct profiler prof;
    theDomain.prof = &prof;
    start_clock( &theDomain ); 
    read_par_file( &theDomain );
  
#if TYPE == 0
    setupDomain_aos( &theDomain );
#endif

    double hash = 0.0;

#if TYPE == 0
    hash = hash_aos(&theDomain);
#endif
    
    printf("HASH: %.16le\n", hash);

    while( theDomain.current_step < theDomain.final_step )
    {
       printf("step: %d\n", theDomain.current_step);
       prof_tick(&prof, PROF_TIMESTEP);
       timestep(&theDomain);
       prof_tock(&prof, PROF_TIMESTEP);
    }

#if TYPE == 0
    hash = hash_aos(&theDomain);
#endif

    printf("HASH: %.16le\n", hash);

    generate_log( &theDomain );
#if TYPE == 0
    freeDomain_aos( &theDomain );
#endif

    return 0;

}
