
#include "header.h"
#include "profiler.h"

void start_clock( struct domain * theDomain )
{
    theDomain->Wallt_init = time(NULL);
    prof_init(theDomain->prof);
    prof_tick(theDomain->prof, PROF_TOT);
}
 
int count_cells( struct domain * theDomain )
{
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;

    int Nc=0;
    int jmin = NgRa;
    int jmax = Nr-NgRb;
    int kmin = NgZa;
    int kmax = Nz-NgZb;

    int j,k;
    for( j=jmin ; j<jmax ; ++j )
        for( k=kmin ; k<kmax ; ++k )
            Nc += Np[j+Nr*k];
   
    return(Nc);
}

void generate_log( struct domain * theDomain )
{
    time_t endtime = time(NULL);
    prof_tock(theDomain->prof, PROF_TOT);
    int seconds = (int) (endtime - theDomain->Wallt_init);

   
    int Nc = count_cells( theDomain );
    int Nt = theDomain->current_step;

    double avgdt = (double)seconds/2./(double)Nc/(double)Nt;

    double prof_secs[NUM_PROF];
    int i;
    for(i=0; i<NUM_PROF; i++)
        prof_secs[i] = ((double)theDomain->prof->elapsed_ticks[i])
                        / CLOCKS_PER_SEC;

    FILE * logfile = fopen("times.log","w");
    fprintf(logfile,"Total time = %d sec\n",seconds);
    fprintf(logfile,"Number of cells = %d\n",Nc);
    fprintf(logfile,"Number of timesteps = %d\n", Nt);
    fprintf(logfile,"Megazones per second = %.2e\n",1./(avgdt*1e6));
    fprintf(logfile,"Megazones per CPU second = %.2e\n",
            1./(avgdt*1e6));
    fprintf(logfile,"Time/zone/step = %.2e microseconds\n",(avgdt*1e6));

    fprintf(logfile, "\n");


    fprintf(logfile, "Profiler total:    %lf sec\n", prof_secs[PROF_TOT]);
    fprintf(logfile, "dt:       %lf sec (%.2lf%%)\n",
            prof_secs[PROF_DT],
            prof_secs[PROF_DT]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile, "timestep: %lf sec (%.2lf%%)\n",
            prof_secs[PROF_TIMESTEP], 
            prof_secs[PROF_TIMESTEP]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "    recon:     %lf sec (%.2lf%% of timestep) (%.2lf%%)\n",
            prof_secs[PROF_RECON], 
            prof_secs[PROF_RECON]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_RECON]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "        recon_r:     %lf sec (%.2lf%% recon, %.2lf%% timestep, %.2lf%%)\n",
            prof_secs[PROF_RECON_R], 
            prof_secs[PROF_RECON_R]*100.0/prof_secs[PROF_RECON],
            prof_secs[PROF_RECON_R]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_RECON_R]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "        recon_p:     %lf sec (%.2lf%% recon, %.2lf%% timestep, %.2lf%%)\n",
            prof_secs[PROF_RECON_P], 
            prof_secs[PROF_RECON_P]*100.0/prof_secs[PROF_RECON],
            prof_secs[PROF_RECON_P]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_RECON_P]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "        recon_z:     %lf sec (%.2lf%% recon, %.2lf%% timestep, %.2lf%%)\n",
            prof_secs[PROF_RECON_Z], 
            prof_secs[PROF_RECON_Z]*100.0/prof_secs[PROF_RECON],
            prof_secs[PROF_RECON_Z]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_RECON_Z]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "    flux:      %lf sec (%.2lf%% of timestep) (%.2lf%%)\n",
            prof_secs[PROF_FLUX], 
            prof_secs[PROF_FLUX]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_FLUX]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "        flux_r:     %lf sec (%.2lf%% flux, %.2lf%% timestep, %.2lf%%)\n",
            prof_secs[PROF_FLUX_R], 
            prof_secs[PROF_FLUX_R]*100.0/prof_secs[PROF_FLUX],
            prof_secs[PROF_FLUX_R]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_FLUX_R]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "        flux_p:     %lf sec (%.2lf%% flux, %.2lf%% timestep, %.2lf%%)\n",
            prof_secs[PROF_FLUX_P], 
            prof_secs[PROF_FLUX_P]*100.0/prof_secs[PROF_FLUX],
            prof_secs[PROF_FLUX_P]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_FLUX_P]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "        flux_z:     %lf sec (%.2lf%% flux, %.2lf%% timestep, %.2lf%%)\n",
            prof_secs[PROF_FLUX_Z], 
            prof_secs[PROF_FLUX_Z]*100.0/prof_secs[PROF_FLUX],
            prof_secs[PROF_FLUX_Z]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_FLUX_Z]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "    ct:        %lf sec (%.2lf%% of timestep) (%.2lf%%)\n",
            prof_secs[PROF_CT], 
            prof_secs[PROF_CT]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_CT]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "    source:    %lf sec (%.2lf%% of timestep) (%.2lf%%)\n",
            prof_secs[PROF_SOURCE], 
            prof_secs[PROF_SOURCE]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_SOURCE]*100.0/prof_secs[PROF_TOT]);
    fprintf(logfile,
            "    cons2prim: %lf sec (%.2lf%% of timestep) (%.2lf%%)\n",
            prof_secs[PROF_C2P], 
            prof_secs[PROF_C2P]*100.0/prof_secs[PROF_TIMESTEP],
            prof_secs[PROF_C2P]*100.0/prof_secs[PROF_TOT]);

    fclose(logfile);
}

void prof_init(struct profiler *prof)
{
    int i;
    for(i=0; i<NUM_PROF; i++)
    {
        prof->ticks[i] = 0;
        prof->elapsed_ticks[i] = 0;
    }
}

void prof_tick(struct profiler *prof, int label)
{
    prof->ticks[label] = clock();
}

void prof_tock(struct profiler *prof, int label)
{
    prof->elapsed_ticks[label] += clock() - prof->ticks[label];
}

