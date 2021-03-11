
#include "header.h"
#include "grid.h"
#include "domain_aos.h"
#include "domain_soa1.h"
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
#endif

    double hash = 0.0;

#if TYPE == 0
    hash = hash_aos(&theDomain, 0);
    printf("HASH0: %.16le\n", hash);
    hash = hash_aos(&theDomain, 1);
    printf("HASH1: %.16le\n", hash);
    hash = hash_aos(&theDomain, 2);
    printf("HASH2: %.16le\n", hash);
    hash = hash_aos(&theDomain, 3);
    printf("HASH3: %.16le\n", hash);
    hash = hash_aos(&theDomain, 4);
    printf("HASH4: %.16le\n", hash);
    hash = hash_aos(&theDomain, 5);
    printf("HASH5: %.16le\n", hash);
    hash = hash_aos(&theDomain, 6);
    printf("HASH6: %.16le\n", hash);
    hash = hash_aos(&theDomain, -1);
    printf("HASH: %.16le\n", hash);
#elif TYPE == 1
    hash = hash_soa1(&theDomain, 0);
    printf("HASH0: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 1);
    printf("HASH1: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 2);
    printf("HASH2: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 3);
    printf("HASH3: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 4);
    printf("HASH4: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 5);
    printf("HASH5: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 6);
    printf("HASH6: %.16le\n", hash);
    hash = hash_soa1(&theDomain, -1);
    printf("HASH: %.16le\n", hash);
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
#if TYPE == 0
    hash = hash_aos(&theDomain, 0);
    printf("HASH0: %.16le\n", hash);
    hash = hash_aos(&theDomain, 1);
    printf("HASH1: %.16le\n", hash);
    hash = hash_aos(&theDomain, 2);
    printf("HASH2: %.16le\n", hash);
    hash = hash_aos(&theDomain, 3);
    printf("HASH3: %.16le\n", hash);
    hash = hash_aos(&theDomain, 4);
    printf("HASH4: %.16le\n", hash);
    hash = hash_aos(&theDomain, 5);
    printf("HASH5: %.16le\n", hash);
    hash = hash_aos(&theDomain, 6);
    printf("HASH6: %.16le\n", hash);
    hash = hash_aos(&theDomain, -1);
    printf("HASH: %.16le\n", hash);
#elif TYPE == 1
    hash = hash_soa1(&theDomain, 0);
    printf("HASH0: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 1);
    printf("HASH1: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 2);
    printf("HASH2: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 3);
    printf("HASH3: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 4);
    printf("HASH4: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 5);
    printf("HASH5: %.16le\n", hash);
    hash = hash_soa1(&theDomain, 6);
    printf("HASH6: %.16le\n", hash);
    hash = hash_soa1(&theDomain, -1);
    printf("HASH: %.16le\n", hash);
#endif

    printf("HASH: %.16le\n", hash);

    generate_log( &theDomain );
#if TYPE == 0
    freeDomain_aos( &theDomain );
#elif TYPE == 1
    freeDomain_soa1( &theDomain );
#endif

    return 0;

}
