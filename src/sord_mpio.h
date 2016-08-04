void read_time_series(char *fname, int nx, int ny, int xi, int yi, 
                        int nt, int nsites, float *buf);
                        
void read_fault_params(char *fname, MPI_Offset off, int nsites, float *buf);


void write_scalar_int(char *fname, MPI_Offset off, int count, int *buf);

/* write out fast component, then time 
ie. mxx(1), myy(1), mzz(1), mxz(1), myz(1), mzz(1), mxx(2) .... mxx(n) */
void write_tensor_time_series(char *fname, MPI_Offset off, int nt, float **buf);

void write_momrate(char* fname, int nst, int nchunks, int rank, int csize, int ion, int *xi, int *yi, int *zi,
     float **xx, float **yy, float **zz, float **xz, float **yz, float **xy);

void write_locations(char *fname, MPI_Offset off, int nsites, int *xloc, int *yloc, int *zloc);
