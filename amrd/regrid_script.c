//=============================================================================
// regrid_script.c --- reading/writing regridding script
//=============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "globals.h"
#include "regrid_script.h"
#include "util.h"

FILE *rgs_stream=0;

//=============================================================================
// Reads the next bit of regridding info from the specified file.
// Returns 1 if OK.
//
// returns:
//
// t    : time of regrid
// L1   : coarsest level overwhich hierarchy can change
// Lf   : finest level of new hierarhcy
// gnum : number of bboxes
// glev : level #'s
// gbbox : bbox's       
//
// caller must free() gnum,glev,gbbox.
//=============================================================================
int rgs_read_next(real *t, int *L1, int *Lf, int *gnum, int **glev, real **gbbox)
{
   static int first=1;
   int i,ret;
   char buffer[1024];
   int ltrace=0;

   if (AMRD_regrid_script!=AMRD_READ_REGRID_SCRIPT) return 0;

   if (first)
   {
      if (!AMRD_regrid_script_name) AMRD_stop("rgs_read_next: error ... no regrid script specified","");
         
      if (!(rgs_stream=fopen(AMRD_regrid_script_name,"r"))) 
         AMRD_stop("rgs_read_next: error ... cannot open file:",AMRD_regrid_script_name);
   }

   if (fscanf(rgs_stream,"%s\n",buffer)<=0) { return 0;}
   if (fscanf(rgs_stream,"%lf,%i,%i,%i\n",t,L1,Lf,gnum)<=0) return 0;

   IFL0 printf("\nrgs_read_next: t=%lf, L1=%i, Lf=%i, gnum=%i\n\n",*t,*L1,*Lf,*gnum);

   if (*gnum<0)
   {
      AMRD_stop("rgs_read_next: error ... gnum < 0","");
   }
   if (*gnum)
   {
       if (!(*glev=(int *)malloc(sizeof(int)*(*gnum)))) AMRD_stop("rgs_read_next: error ... out of memory","");
       if (!(*gbbox=(real *)malloc(2*AMRD_dim*sizeof(real)*(*gnum)))) AMRD_stop("rgs_read_next: error ... out of memory","");
   }
   else
   {
      *glev=0;
      *gbbox=0;
   }

   for (i=0; i<*gnum; i++)
   {
      switch(AMRD_dim)
      {
         case 1: ret=fscanf(rgs_stream,"%i,%lf,%lf\n",&((*glev)[i]),&((*gbbox)[2*i]),&((*gbbox)[2*i+1])); break;
         case 2: ret=fscanf(rgs_stream,"%i,%lf,%lf,%lf,%lf\n",&((*glev)[i]),
                 &((*gbbox)[4*i]),&((*gbbox)[4*i+1]),&((*gbbox)[4*i+2]),&((*gbbox)[4*i+3])); break;
         case 3: ret=fscanf(rgs_stream,"%i,%lf,%lf,%lf,%lf,%lf,%lf\n",&((*glev)[i]),
                 &((*gbbox)[6*i]),&((*gbbox)[6*i+1]),&((*gbbox)[6*i+2]),&((*gbbox)[6*i+3]),&((*gbbox)[6*i+4]),
                 &((*gbbox)[6*i+5])); break;
         default: AMRD_stop("rgs_read_next: specified dim not suppported\n","");
      }
      if (ret<=0) return 0;
   }

   first=0;
   return 1;
}

//=============================================================================
// The 'inverse' of rgs_read_next()
//
// NOTE: if this was a restart, then open file in mode 'append'
//=============================================================================
int rgs_write_next(real t, int L1, int Lf, int gnum, int *glev, real *gbbox)
{
   static int first=1;
   int i,ltrace=0;

   if (AMRD_regrid_script!=AMRD_WRITE_REGRID_SCRIPT) return 0;

   if (my_rank!=0) return 1;

   IFL printf("\nrgs_write_next: t=%lf, L1=%i, Lf=%i, gnum=%i\n\n",t,L1,Lf,gnum);

   if (first)
   {
      if (!AMRD_regrid_script_name) AMRD_stop("rgs_write_next: error ... no regrid script specified\n","");
        
      if (AMRD_cp_restart)
      {
         if (!my_rank) printf("NOTE: appending to regrid script file, as restart flag is on\n"); 
         if (!(rgs_stream=fopen(AMRD_regrid_script_name,"a"))) 
            AMRD_stop("rgs_write_next: error ... cannot open file:",AMRD_regrid_script_name);
      }
      else
      {
         if (!(rgs_stream=fopen(AMRD_regrid_script_name,"w"))) 
            AMRD_stop("rgs_write_next: error ... cannot open file:",AMRD_regrid_script_name);
      }
   }

   fprintf(rgs_stream,"Regrid:t=%lf,L1=%i,Lf=%i,gnum=%i\n",t,L1,Lf,gnum);
           fprintf(rgs_stream,"%lf,%i,%i,%i\n",t,L1,Lf,gnum);
   for (i=0; i<gnum; i++)
   {
      switch(AMRD_dim)
      {
         case 1: fprintf(rgs_stream,"%i,%21.15lf,%21.15lf\n",glev[i],gbbox[2*i],gbbox[2*i+1]); break;
         case 2: fprintf(rgs_stream,"%i,%21.15lf,%21.15lf,%21.15lf,%21.15lf\n",glev[i],
                 gbbox[4*i],gbbox[4*i+1],gbbox[4*i+2],gbbox[4*i+3]); break;
         case 3: fprintf(rgs_stream,"%i,%21.15lf,%21.15lf,%21.15lf,%21.15lf,%21.15lf,%21.15lf\n",glev[i],
                 gbbox[6*i],gbbox[6*i+1],gbbox[6*i+2],gbbox[6*i+3],gbbox[6*i+4],gbbox[6*i+5]); break;
         default: AMRD_stop("rgs_write_next: specified dim not suppported\n","");
      }
   }

   fflush(rgs_stream);

   first=0;
   return 1;
}
