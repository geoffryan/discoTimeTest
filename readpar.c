enum{VAR_INT,VAR_DOUB,VAR_STR};

#include "header.h"

int readvar( char * filename , char * varname , int vartype , void * ptr ){

   FILE * inFile = fopen( filename , "r" );
   char s[512];
   char nm[512]="";
   char s1[512];
   int found = 0;
   
   while( (fgets(s,512,inFile) != NULL) && found==0 ){
      sscanf(s,"%s ",nm);
      if( strcmp(nm,varname)==0 ){
         strcpy(s1,s);
         found=1;
      }
   }
   
   fclose( inFile );
   if( found==0 ) return(1);

   char * s2 = s1+strlen(nm)+strspn(s1+strlen(nm),"\t :=>_");

   double temp;
   char stringval[256];

   sscanf(s2,"%lf",&temp);
   sscanf(s2,"%256s",stringval);

   if( vartype == VAR_INT ){
      *((int *)   ptr) = (int)temp;
   }else if( vartype == VAR_DOUB ){
      *((double *)ptr) = (double)temp;
   }else{
      strcpy( ptr , stringval );
   }

   return(0);
}

int read_par_file(struct domain * theDomain)
{
    struct param_list * theList = &( theDomain->theParList );

    char pfile[] = "in.par";

    int err = 0;  

    err += readvar(pfile, "Num_R",   VAR_INT, &(theList->Num_R));
    err += readvar(pfile, "Num_Z",   VAR_INT, &(theList->Num_Z));
    err += readvar(pfile, "nsteps",  VAR_INT, &(theList->nsteps));
    err += readvar(pfile, "R_Min",   VAR_DOUB, &(theList->rmin));
    err += readvar(pfile, "R_Max",   VAR_DOUB, &(theList->rmax));
    err += readvar(pfile, "Z_Min",   VAR_DOUB, &(theList->zmin));
    err += readvar(pfile, "Z_Max",   VAR_DOUB, &(theList->zmax));
    err += readvar(pfile, "Phi_Max", VAR_DOUB, &(theList->phimax));
    err += readvar(pfile, "dt", VAR_DOUB, &(theList->dt));

    theList->phimax *= 2*M_PI;

    if( err > 0 )
    {
        printf("Read Failed\n");
        return err;
    }

    return 0;
}


