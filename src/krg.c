/* This program creates a kinematic source from results of an AWP dynamic rupture
   simulation, both for a vertical non-planar fault.

    William Savran, wsavran@ucsd.edu

   MPI-IO is used both for reading the source time function and writing
   the moment rate file.
*/

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include "sord_mpio.h"
#include "utils.h"
#include "xapiir.h"
#include "spline.h"
#include "params.h"

/* computed scalar moment given moment rate tensor */
float calculate_moment(int nst, float dt, float *xx, float *yy, float *zz, float *xz, float *yz, float *xy);

/* calculate moment rate tensor given subfault parameters */
void calculate_moment_rate(float *time, int nt, float *xx, float *yy, float *zz, float *xz, float *yz, float *xy, 
                                    float moment, float strike, float dip, float rake, float psv, float trup);

int main (int argc, char*argv[]) {
    /* modify these parameters */
    param p;
    read_settings(&p, "./in/params.txt");
    int s0;
    MPI_Offset off;
    int xi, yi;
    int rank, nprocs, csize;
    int k,l;
    float moment, global_moment;
    float *psv_buf, *trup_buf;
    float *strike_buf, *dip_buf, *rake_buf;
    float *coords_buf;
    float *moment_buf;
    int debug = 0;
    int master = 0;
    float **xx_buf, **yy_buf, **zz_buf, **xz_buf, **yz_buf, **xy_buf;
    int *xil_buf, *yil_buf, *zil_buf;
    float *time_buf;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (debug == 1) {
        // output some parameters
        fprintf(stderr, "rank: %i\n", rank); 
        fprintf(stderr, "s0: %i\n", s0);
        fprintf(stderr, "off: %lli\n", off);
        fprintf(stderr, "yi: %i\n", yi);
        fprintf(stderr, "xi: %i\n", xi);
        fprintf(stderr, "nx: %i\n", p.nx);
        fprintf(stderr, "nz: %i\n", p.nz);
        fprintf(stderr, "dt: %f\n", p.dt);
        fprintf(stderr, "dx: %f\n", p.dx);
        fprintf(stderr, "nt: %i\n", p.nt);
        fprintf(stderr, "psv file: %s\n", p.psv_file);
        fprintf(stderr, "trup file: %s\n", p.trup_file);
        fprintf(stderr, "strike file: %s\n", p.strike_file);
        fprintf(stderr, "dip file: %s\n", p.dip_file);
        fprintf(stderr, "rake file: %s\n", p.rake_file);
        fprintf(stderr, "moment file: %s\n", p.moment_file);
        fprintf(stderr, "momentrate file: %s\n\n", p.momentrate_file);
    }

    /* number of points to process per read operation */
    if (((p.nx*p.nz) % (nprocs * p.nchunks)) != 0) {
        fprintf(stdout, "number of points not divisible by number of cpus * buffer size\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(1);
    } 
    
    /* number of time-series read per mpi-io call */
    csize = p.nx*p.nz / nprocs / p.nchunks;
    
    time_buf = arange(0.0, p.rt, p.dt, &(p.nt));
        /* allocate arrays */
    psv_buf = (float*)calloc(csize, sizeof(float));
    trup_buf = (float*)calloc(csize, sizeof(float));
    strike_buf = (float*)calloc(csize, sizeof(float));
    dip_buf = (float*)calloc(csize, sizeof(float));
    rake_buf = (float*)calloc(csize, sizeof(float));
    coords_buf = (float*)calloc(csize, sizeof(float));
    moment_buf = (float*)calloc(csize, sizeof(float));

    time_buf = (float*)calloc(p.nt, sizeof(float));

    xil_buf = (int*)calloc(csize, sizeof(int));
    yil_buf = (int*)calloc(csize, sizeof(int));
    zil_buf = (int*)calloc(csize, sizeof(int));


    /* 2D arrays for time-series */
    //xx_buf = (float**)calloc(csize, sizeof(float*));
    //yy_buf = (float**)calloc(csize, sizeof(float*));
    //zz_buf = (float**)calloc(csize, sizeof(float*));
    //xz_buf = (float**)calloc(csize, sizeof(float*));
    //yz_buf = (float**)calloc(csize, sizeof(float*));
    //xy_buf = (float**)calloc(csize, sizeof(float*));
    //
    //for (l=0; l<csize; l++) {
    //    xx_buf[l] = (float*)calloc(p.nt, sizeof(float));
    //    yy_buf[l] = (float*)calloc(p.nt, sizeof(float));
    //    zz_buf[l] = (float*)calloc(p.nt, sizeof(float));
    //    xz_buf[l] = (float*)calloc(p.nt, sizeof(float));
    //    yz_buf[l] = (float*)calloc(p.nt, sizeof(float));
    //    xy_buf[l] = (float*)calloc(p.nt, sizeof(float));
    //}

        /* loop over subfault block */
        MPI_Barrier(MPI_COMM_WORLD);
        moment = 0;
        for (l = 0; l < p.nchunks; l++) {
            
            // calculate offsets
            s0 = rank*p.nchunks*csize + l*csize;
            off = (MPI_Offset) s0 * sizeof(float);
            yi = s0 / p.nx;
            xi = s0 % p.nx;   
            if (debug == 0) {
                fprintf(stderr, "s0: %d\n", s0);
                fprintf(stderr, "csize: %d\n", csize);
                fprintf(stderr, "off: %lli\n", off);
            }

            // read files
            read_fault_params(p.moment_file, off, csize, moment_buf);
            read_fault_params(p.psv_file, off, csize, psv_buf);
            read_fault_params(p.trup_file, off, csize, trup_buf);
            read_fault_params(p.strike_file, off, csize, strike_buf);
            read_fault_params(p.dip_file, off, csize, dip_buf);
            read_fault_params(p.rake_file, off, csize, rake_buf);
            read_fault_params(p.coord_file, off, csize, coords_buf);

            if (rank == 0) {
                for (k=0; k<csize; k++) {
                    fprintf(stdout, "%f\n", moment_buf[k]);
                }

            }

            /* looping over each subfault */
           // for (k=0; k<csize; k++) {
           //     calculate_moment_rate(time_buf, p.nt, xx_buf[k], yy_buf[k], zz_buf[k], xz_buf[k], yz_buf[k], xy_buf[k],
           //                             moment_buf[k], strike_buf[k], dip_buf[k], rake_buf[k], psv_buf[k], trup_buf[k]);

           //     /* determine subfault location assuming fault is striking along x component */
           //     xil_buf[k] = p.x_start + (s0+k) % p.nx; 
           //     yil_buf[k] = p.y_start + rint((coords_buf[k] - p.mean_faultn_coord) / p.dx);
           //     zil_buf[k] = p.z_start + (s0+k) / p.nx;  

           //     moment += calculate_moment(p.nt, p.dt, xx_buf[k], yy_buf[k], zz_buf[k], xz_buf[k], yz_buf[k], xy_buf[k]);
           // }
           // 
           // if (rank==master) fprintf(stderr, "writing momentrate file.");
//         // write_momrate(p.momentrate_file, p.nt, p.nchunks, rank, csize, l, xil_buf, yil_buf, zil_buf, xx_buf, yy_buf, zz_buf, xz_buf, yz_buf, xy_buf);
        
    } /* end main loop */

    /* compute moment from all processes */
    //MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Reduce(&moment, &global_moment, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
    //if (rank == master) {
    //    fprintf(stderr, "total moment: %e\n", global_moment);
    //    fprintf(stderr, "mw: %e\n", (2.0/3.0)*(log10(global_moment)-9.1));
    //}
    /* free buffers */
    
    /* finalize mpi */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

void calculate_moment_rate(float *time, int nt, float *xx, float *yy, float *zz, float *xz, float *yz, float *xy, 
                                    float moment, float strike, float dip, float rake, float psv, float trup) {

    float momentrate;
    float tpeak;
    int i;

    tpeak = moment / exp(psv);
    for (i=0; i<nt; i++) {
        if (time[i] > trup) {
            momentrate = moment / tpeak * (time[i] - trup) / tpeak * exp(-((time[i] - trup) / tpeak));
        } else {
            momentrate = 0.0;
        }
        // convert to magnitude from strike dip and rake
        xx[i] = momentrate * (sin(dip)*cos(rake)*sin(2*strike) - sin(2*dip)*sin(rake)*cos(strike)*cos(strike));
        yy[i] = -momentrate * (sin(dip)*cos(rake)*sin(2*strike) + sin(2*dip)*sin(rake)*sin(strike)*sin(strike));
        zz[i] = momentrate * (sin(2*dip)*sin(rake));
        xz[i] = momentrate * (cos(dip)*cos(rake)*sin(strike) - cos(2*dip)*sin(rake)*cos(strike));
        yz[i] = momentrate * (cos(dip)*cos(rake)*cos(strike)+cos(2*dip)*sin(rake)*sin(strike));
        xy[i] = momentrate * (sin(dip)*cos(rake)*cos(2*strike) + 0.5*sin(2*dip)*sin(rake)*sin(2*strike));
    }
}
     
float calculate_moment(int nst, float dt, float *xx, float *yy, float *zz, float *xz, float *yz, float *xy) {
    float mxx, myy, mzz, mxz, myz, mxy;
    float mij_squared;
    float m0;
    float prefactor;

    /* integrate moment-rate to get moment */
    mxx = sum(xx, nst)*dt;
    myy = sum(yy, nst)*dt;
    mzz = sum(zz, nst)*dt;
    mxz = sum(xz, nst)*dt;
    myz = sum(yz, nst)*dt;
    mxy = sum(xy, nst)*dt;

    /* square moment rate tensor */
    prefactor = 1.0f / sqrt(2);
    mij_squared = mxx*mxx + myy*myy + mzz*mzz + 2*mxz*mxz + 2*myz*myz + 2*mxy*mxy;
    m0 = prefactor * sqrt(mij_squared);
    return m0;
}
