#include "params.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void read_settings(param *P, char *paramsfile) {
  /* open the two configuration files*/
   FILE *fid = fopen(paramsfile, "r");
   if (fid==NULL) {
           fprintf(stdout, "%s\n", paramsfile);
	   perror("could not open input file");
	   exit(1);
   }

   /* read the parameters from IN3D */
   char *stat2;
   char param[15], value[100];
   char newline[200];
 
   stat2=(char*) calloc(2, sizeof(char));
   while (stat2 != NULL){
     stat2=fgets(newline, 200, fid);
     sscanf(newline, "%s %s", param, value);
     if (strcmp(param, "nchunks")==0) sscanf(value, "%d", &P->nchunks);
     if (strcmp(param, "nx")==0) sscanf(value, "%d", &P->nx);
     if (strcmp(param, "nz")==0) sscanf(value, "%d", &P->nz);
     if (strcmp(param, "dx")==0) sscanf(value, "%f", &P->dx);
     if (strcmp(param, "dt")==0) sscanf(value, "%f", &P->dt);
     if (strcmp(param, "rt")==0) sscanf(value, "%f", &P->rt);
     if (strcmp(param, "nt")==0) sscanf(value, "%d", &P->nt);
     if (strcmp(param, "x_start")==0) sscanf(value, "%d", &P->x_start);
     if (strcmp(param, "y_start")==0) sscanf(value, "%d", &P->y_start);
     if (strcmp(param, "z_start")==0) sscanf(value, "%d", &P->z_start);
     if (strcmp(param, "mean_faultn_coord")==0) sscanf(value, "%f", &P->mean_faultn_coord);
     if (strcmp(param, "nt")==0) sscanf(value, "%d", &P->nt);
     if (strcmp(param, "psv_file")==0)  strcpy(P->psv_file, value);
     if (strcmp(param, "trup_file")==0)  strcpy(P->trup_file, value);
     if (strcmp(param, "strike_file")==0) strcpy(P->strike_file, value);
     if (strcmp(param, "dip_file")==0) strcpy(P->dip_file, value);
     if (strcmp(param, "rake_file")==0) strcpy(P->rake_file, value);
     if (strcmp(param, "moment_file")==0) strcpy(P->moment_file, value);
     if (strcmp(param, "momentrate_file")==0) strcpy(P->momentrate_file, value);
     if (strcmp(param, "coord_file")==0) strcpy(P->coord_file, value);
     if (strcmp(param, "iord")==0) sscanf(value, "%d", &P->iord);
     if (strcmp(param, "npas")==0) sscanf(value, "%d", &P->npas);
     if (strcmp(param, "trbndw")==0) sscanf(value, "%f", &P->trbndw);
     if (strcmp(param, "a")==0) sscanf(value, "%f", &P->a);
     if (strcmp(param, "aproto")==0) strcpy(P->aproto, value);
     if (strcmp(param, "ftype")==0) strcpy(P->ftype, value);
     if (strcmp(param, "hp")==0) sscanf(value, "%f", &P->hp);
     if (strcmp(param, "lp")==0) sscanf(value, "%f", &P->lp);
   }
   fclose(fid);
}

