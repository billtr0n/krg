#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

void error_check(int ierr, char *message) {
   char errmsg[500];
   int errlen;
   if (ierr != MPI_SUCCESS) {
      fprintf(stderr, "%d: Error in %s\n", ierr, message);
      MPI_Error_string(ierr, errmsg, &errlen);
      fprintf(stderr, "%s", errmsg);
      exit(1);
   }
}

void read_time_series(char *fname, int nx, int ny, int xi, int yi, 
                        int nt, int nsites, float *buf) {
    MPI_File fh;
    MPI_Datatype filetype;
    MPI_Offset disp;
    int ierr;
    
    ierr=MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    error_check(ierr, "MPI_File_open()");
    
    ierr=MPI_Type_vector(nt, nsites, nx*ny, MPI_FLOAT, &filetype);
    error_check(ierr, "MPI_Type_vector()");
    
    ierr=MPI_Type_commit(&filetype);
    error_check(ierr, "MPI_Type_commit()");
    
    disp = (MPI_Offset) sizeof(float) * (nx*yi + xi);
    ierr = MPI_File_set_view(fh, disp, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    error_check(ierr, "MPI_FILE_read_all()");
    
    ierr=MPI_File_read_all(fh, buf, nt*nsites, MPI_FLOAT, MPI_STATUS_IGNORE);
    error_check(ierr, "MPI_File_read_all()");
    
    MPI_File_close(&fh);
    
    ierr=MPI_Type_free(&filetype);
    error_check(ierr, "MPI_Type_free()");
}

void read_fault_params(char *fname, MPI_Offset off, int nsites, float *buf) {
    MPI_File fh;
    int ierr;

    ierr=MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    error_check(ierr, "MPI_File_open()");    
    
    ierr=MPI_File_read_at_all(fh, off, buf, nsites, MPI_FLOAT, MPI_STATUS_IGNORE);
    error_check(ierr, "MPI_File_read_all()");  
    
    MPI_File_close(&fh);
}

void write_locations(char *fname, MPI_Offset off, int nsites, int *xloc, int *yloc, int *zloc) {
    MPI_File fh;
    int ierr;
    int *buf;
    int i;

    buf = (int*)malloc(sizeof(int)*nsites*3);

    ierr=MPI_File_open(MPI_COMM_WORLD, "tpsrc_tmp", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    error_check(ierr, "MPI_File_open()");

    // fill buffer with x,y,z locations 
    for (i=0; i<nsites; i++) {
        buf[3*i] = xloc[i];
        buf[3*i+1] = yloc[i];
        buf[3*i+2] = zloc[i];
    }

    ierr = MPI_File_write_at_all(fh, 3*off, buf, 3*nsites, MPI_INT, MPI_STATUS_IGNORE);
    error_check(ierr, "MPI_File_write_at_all()");

    MPI_File_close(&fh);
}
/* non-collective i/o ... DEPRECATED */
void write_scalar_int(char *fname, MPI_Offset off, int count, int *buf) {
    MPI_File fh;
    int ierr;
    
    ierr=MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    error_check(ierr, "MPI_File_open()");    

    ierr=MPI_File_write_at(fh, off, buf, count, MPI_INT, MPI_STATUS_IGNORE);
    error_check(ierr, "MPI_File_write_at_all_1()"); 
    
    MPI_File_close(&fh);    
}


/* non-collective i/o ... DEPRECATED */
void write_tensor_time_series(char *fname, MPI_Offset off, int nt, float **buf) {
    MPI_File fh;
    int ierr;
    int l, n;
    float *out_buf;
    int count;
    
    count = nt * 6;
    out_buf = (float*)malloc(nt*6*sizeof(float)); /* 6 for symmetric tensor */
    
    /* map buf to 1d array */
    for (l=0; l<6; l++) {
        for (n=0; n<nt; n++) {
            out_buf[n*6+l] = buf[l][n];
        }
    }

    ierr=MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    error_check(ierr, "MPI_File_open()");    

    ierr=MPI_File_write_at(fh, off, out_buf, count, MPI_FLOAT, MPI_STATUS_IGNORE);
    error_check(ierr, "MPI_File_write_at()"); 
    
    MPI_File_close(&fh);    
    free(out_buf);
} 

void write_momrate(char *fname, int nst, int nchunks, int rank, int csize, int ion, int *xi, int *yi, int *zi,
     float **xx, float **yy, float **zz, float **xz, float **yz, float **xy){

   MPI_Offset offset;
   int *buf1, *blen1, *blen2;
   int ierr;
   int p, n, bpos;
   float *buf2;

   MPI_Aint *map1,*map2;
   MPI_File fh;
   MPI_Datatype filetype1, filetype2;

   MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

   buf1 = (int*) calloc(3*csize, sizeof(int));
   map1 = (MPI_Aint*) calloc(csize, sizeof(MPI_Aint));
   blen1 = (int*) calloc(csize, sizeof(int));

   for (p=0; p<csize; p++) {
      // map1[p]=(MPI_Aint) (3*sizeof(int) + nst*6*sizeof(float))  * (MPI_Aint) + p;
      map1[p]=(MPI_Aint) (3*sizeof(int) + nst*6*sizeof(float))  * (MPI_Aint) p;
      //if (rank==0) fprintf(stdout, "map[%d]=%ld\n", p, map1[p]);
      blen1[p]=3;
      buf1[p*3]=xi[p];
      buf1[p*3+1]=yi[p];
      buf1[p*3+2]=zi[p];
   }
   
   offset = (MPI_Offset) (3*sizeof(int) + nst*6*sizeof(float)) 
      * (MPI_Offset) (csize*rank*nchunks+ion*csize);
   if (rank==0) fprintf(stdout, "offset = %lld\n", offset);
   if (rank==0) fprintf(stdout, "ion = %d, csize=%d\n", ion, csize);

   ierr=MPI_Type_create_hindexed(csize, blen1, map1, MPI_INT, &filetype1);
   error_check(ierr, "MPI_Type_create_hindexed()");

   ierr=MPI_Type_commit(&filetype1);
   error_check(ierr, "MPI_Type_commit()");

   ierr=MPI_File_set_view(fh, offset, MPI_INT, filetype1, "native", MPI_INFO_NULL);
   error_check(ierr, "MPI_File_set_view()");
   
   ierr=MPI_File_write_all(fh, buf1, csize*3, MPI_INT, MPI_STATUS_IGNORE);
   error_check(ierr, "MPI_File_write_all()");

   ierr=MPI_Type_free(&filetype1);
   error_check(ierr, "MPI_Type_free()");

   free(buf1);
   free(map1);
   free(blen1);

   buf2 = (float*) calloc(csize*nst*6, sizeof(float));
   map2 = (MPI_Aint*) calloc(csize, sizeof(MPI_Aint));
   blen2 = (int*) calloc(csize, sizeof(int));

   for (p=0; p<csize; p++){
      map2[p]=(MPI_Aint) (3*sizeof(int) + nst*6*sizeof(float)) * (MPI_Aint) p + 3*sizeof(int);
      blen2[p]=nst*6;

      for (n=0; n<nst; n++) {
         bpos=p*nst*6 + n*6;
         buf2[bpos] = xx[p][n];
         buf2[bpos+1] = yy[p][n];
         buf2[bpos+2] = zz[p][n];
         buf2[bpos+3] = xz[p][n];
         buf2[bpos+4] = yz[p][n];
         buf2[bpos+5] = xy[p][n];
      }
   }

   ierr=MPI_Type_create_hindexed(csize, blen2, map2, MPI_FLOAT, &filetype2);
   error_check(ierr, "MPI_Type_create_hindexed()");

   ierr=MPI_Type_commit(&filetype2);
   error_check(ierr, "MPI_Type_commit()");

   ierr=MPI_File_set_view(fh, offset, MPI_FLOAT, filetype2, "native", MPI_INFO_NULL);
   error_check(ierr, "MPI_File_set_view()");

   ierr=MPI_File_write_all(fh, buf2, csize*6*nst, MPI_FLOAT, MPI_STATUS_IGNORE);
   error_check(ierr, "MPI_File_write_at_all()");

   /* test 
   MPI_Barrier(MPI_COMM_WORLD); */

   ierr=MPI_Type_free(&filetype2);
   error_check(ierr, "MPI_Type_free()");

   ierr=MPI_File_close(&fh);
   error_check(ierr, "MPI_File_close()");

   free(buf2);
   free(map2);
   free(blen2);
}

/* test driver for mpi-i /o subroutines */
/*
int main() {
    return 0;
}
*/
