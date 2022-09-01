
#include "header.h"
#include "grid.h"
#include "domain_aos.h"
#include "domain_soa1.h"
#include "domain_soa2.h"
#include "substep_aos.h"
#include "substep_soa1.h"
#include "substep_soa2.h"
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
    setupGrid(&theDomain);
  
#if TYPE == 0
    setupDomain_aos( &theDomain );
#elif TYPE == 1
    setupDomain_soa1( &theDomain );
#elif TYPE == 2
    setupDomain_soa2( &theDomain );
#endif

    double hash = 0.0;

#if TYPE == 0
    hash = hash_aos(&theDomain, 0);
    printf("HASH0:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 1);
    printf("HASH1:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 2);
    printf("HASH2:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 3);
    printf("HASH3:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 4);
    printf("HASH4:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 5);
    printf("HASH5:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 6);
    printf("HASH6:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, -1);
    printf("HASH:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
#elif TYPE == 1
    hash = hash_soa1(&theDomain, 0);
    printf("HASH0:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 1);
    printf("HASH1:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 2);
    printf("HASH2:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 3);
    printf("HASH3:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 4);
    printf("HASH4:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 5);
    printf("HASH5:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 6);
    printf("HASH6:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, -1);
    printf("HASH:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
#elif TYPE == 2
    hash = hash_soa2(&theDomain, 0);
    printf("HASH0:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 1);
    printf("HASH1:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 2);
    printf("HASH2:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 3);
    printf("HASH3:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 4);
    printf("HASH4:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 5);
    printf("HASH5:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 6);
    printf("HASH6:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, -1);
    printf("HASH:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
#endif
    
    while( theDomain.current_step < theDomain.final_step )
    {
       printf("step: %d\n", theDomain.current_step);
       prof_tick(&prof, PROF_TIMESTEP);
       timestep(&theDomain);
       prof_tock(&prof, PROF_TIMESTEP);
    }
/*
#if TYPE == 0
    hash = hash_aos(&theDomain, -1);
#elif TYPE == 1
    hash = hash_soa1(&theDomain, -1);
#endif
*/
    generate_log(&theDomain);

    dump(&theDomain, "fin");

#if TYPE == 0
    hash = hash_aos(&theDomain, 0);
    printf("HASH0:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 1);
    printf("HASH1:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 2);
    printf("HASH2:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 3);
    printf("HASH3:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 4);
    printf("HASH4:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 5);
    printf("HASH5:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, 6);
    printf("HASH6:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_aos(&theDomain, -1);
    printf("HASH:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
#elif TYPE == 1
    hash = hash_soa1(&theDomain, 0);
    printf("HASH0:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 1);
    printf("HASH1:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 2);
    printf("HASH2:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 3);
    printf("HASH3:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 4);
    printf("HASH4:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 5);
    printf("HASH5:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, 6);
    printf("HASH6:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa1(&theDomain, -1);
    printf("HASH:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
#elif TYPE == 2
    hash = hash_soa2(&theDomain, 0);
    printf("HASH0:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 1);
    printf("HASH1:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 2);
    printf("HASH2:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 3);
    printf("HASH3:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 4);
    printf("HASH4:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 5);
    printf("HASH5:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, 6);
    printf("HASH6:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
    hash = hash_soa2(&theDomain, -1);
    printf("HASH:  %lx  %.16le\n", *((unsigned long *) &hash), hash);
#endif

    printf("HASH:  %lx  %.16le\n", *((unsigned long *) &hash), hash);


#if TYPE == 0
    //dump_grid_aos(&theDomain);
    freeDomain_aos(&theDomain);
#elif TYPE == 1
    //dump_grid_soa1(&theDomain);
    freeDomain_soa1(&theDomain);
#elif TYPE == 2
    //dump_grid_soa2(&theDomain);
    freeDomain_soa2(&theDomain);
#endif

    return 0;

}
