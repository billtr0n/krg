/* This program creates a kinematic source from results of an AWP dynamic rupture
   simulation, both for a vertical non-planar fault.  This code was largely based on
   code provided by Daneil Roten <droten at sdsu dot mail dot com>. 

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
#include "stf.h"

/* compute tr based on multiple linear regression */
float compute_tr(float slip, float treff);

/* computed scalar moment given moment rate tensor */
float calculate_moment(int nst, float dt, float *xx, float *yy, float *zz, float *xz, float *yz, float *xy);

/* calculate moment rate tensor given subfault parameters */
void calculate_moment_rate(float *time, int nt, float dx, float *xx, float *yy, float *zz, float *xz, float *yz, float *xy, 
                                    float slip, float psv, float trup, float strike, float dip, float rake, 
                                    float vs, float rho, int rank, float dc, float median_ts, float truptot, float tp_psv_coef);

int main (int argc, char*argv[]) {
    /* modify these parameters */
    param p;
    int s0;
    MPI_Offset off;
    int xi, yi;
    int rank, nprocs, csize;
    int k,l;
    float moment, global_moment;
    float *psv_buf, *trup_buf;
    float *strike_buf, *dip_buf, *rake_buf;
    float *coords_buf;
    float *slip_buf;
    int debug = 0;
    int master = 0;
    float **xx_buf, **yy_buf, **zz_buf, **xz_buf, **yz_buf, **xy_buf;
    float *vs_buf, *rho_buf;
    int *xil_buf, *yil_buf, *zil_buf;
    int nt;
    float *time_buf;
    int source_sim_dx_rat;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* read parameters */
    read_settings(&p, "./in/params.txt");
    source_sim_dx_rat = p.source_dx / p.sim_dx;
    if (source_sim_dx_rat < 1) {
        fprintf(stderr, "source dx must be greater than or equal to sim dx.\n");
        exit(1);
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
    
    /* buffers necessary to create moment-rates */
    time_buf = arange(0.0, p.rt, p.dt, &nt);
    psv_buf = (float*)calloc(csize, sizeof(float));
    vs_buf = (float*)calloc(csize, sizeof(float));
    rho_buf = (float*)calloc(csize, sizeof(float));
    trup_buf = (float*)calloc(csize, sizeof(float));
    strike_buf = (float*)calloc(csize, sizeof(float));
    dip_buf = (float*)calloc(csize, sizeof(float));
    rake_buf = (float*)calloc(csize, sizeof(float));
    coords_buf = (float*)calloc(csize, sizeof(float));
    slip_buf = (float*)calloc(csize, sizeof(float));

    /* subfault locations */
    xil_buf = (int*)calloc(csize, sizeof(int));
    yil_buf = (int*)calloc(csize, sizeof(int));
    zil_buf = (int*)calloc(csize, sizeof(int));

    /* 2D arrays for time-series */
    xx_buf = (float**)calloc(csize, sizeof(float*));
    yy_buf = (float**)calloc(csize, sizeof(float*));
    zz_buf = (float**)calloc(csize, sizeof(float*));
    xz_buf = (float**)calloc(csize, sizeof(float*));
    yz_buf = (float**)calloc(csize, sizeof(float*));
    xy_buf = (float**)calloc(csize, sizeof(float*));
    for (l=0; l<csize; l++) {
        xx_buf[l] = (float*)calloc(nt, sizeof(float));
        yy_buf[l] = (float*)calloc(nt, sizeof(float));
        zz_buf[l] = (float*)calloc(nt, sizeof(float));
        xz_buf[l] = (float*)calloc(nt, sizeof(float));
        yz_buf[l] = (float*)calloc(nt, sizeof(float));
        xy_buf[l] = (float*)calloc(nt, sizeof(float));
    }

    /* output parameters before computation */
    if (rank == master) {
        fprintf(stderr, "nx: %i\n", p.nx);
        fprintf(stderr, "nz: %i\n", p.nz);
        fprintf(stderr, "dt: %f\n", p.dt);
        fprintf(stderr, "sim_dx: %f\n", p.sim_dx);
        fprintf(stderr, "source_dx: %f\n", p.source_dx);
        fprintf(stderr, "rt: %f\n", p.rt);
        fprintf(stderr, "nt: %i\n", nt);
        fprintf(stderr, "proj: %i\n", p.proj);
        fprintf(stderr, "psv file: %s\n", p.psv_file);
        fprintf(stderr, "trup file: %s\n", p.trup_file);
        fprintf(stderr, "strike file: %s\n", p.strike_file);
        fprintf(stderr, "dip file: %s\n", p.dip_file);
        fprintf(stderr, "rake file: %s\n", p.rake_file);
        fprintf(stderr, "slip file: %s\n", p.slip_file);
        fprintf(stderr, "mean fault coord: %f\n", p.faultn_coord);
        fprintf(stderr, "momentrate file: %s\n\n", p.momentrate_file);
        fprintf(stderr, "iord: %d\n", p.iord);
        fprintf(stderr, "npas: %d\n", p.npas);
        fprintf(stderr, "trbndw: %f\n", p.trbndw);
        fprintf(stderr, "a: %f\n", p.a);
        fprintf(stderr, "aproto: %s\n", p.aproto);
        fprintf(stderr, "ftype: %s\n", p.ftype);
        fprintf(stderr, "hp: %f\n", p.hp);
        fprintf(stderr, "lp: %f\n", p.lp);
        fprintf(stderr, "dc: %f\n", p.dc);
        fprintf(stderr, "median_ts: %f\n", p.median_ts);
        fprintf(stderr, "truptot: %f\n", p.truptot);
        fprintf(stderr, "tp_psv_coef: %f\n", p.tp_psv_coef);
    }
    
    /* loop over subfault block nchunks */
    MPI_Barrier(MPI_COMM_WORLD);
    moment = 0;
    for (l = 0; l < p.nchunks; l++) {

        /* calculate offsets */
        s0 = rank*p.nchunks*csize + l*csize;
        off = (MPI_Offset) s0 * sizeof(float);
        yi = s0 / p.nx;
        xi = s0 % p.nx;   
        if (debug > 0) {
            fprintf(stderr, "rank: %d\n", rank);
            fprintf(stderr, "s0: %d\n", s0);
            fprintf(stderr, "csize: %d\n", csize);
            fprintf(stderr, "off: %lli\n", off);
            fprintf(stderr, "xi: %i\n", xi);
            fprintf(stderr, "yi: %i\n", yi);
        }

        /* read files */
        read_fault_params(p.slip_file, off, csize, slip_buf);
        read_fault_params(p.psv_file, off, csize, psv_buf);
        read_fault_params(p.trup_file, off, csize, trup_buf);
        read_fault_params(p.strike_file, off, csize, strike_buf);
        read_fault_params(p.dip_file, off, csize, dip_buf);
        read_fault_params(p.rake_file, off, csize, rake_buf);
        if (p.proj == 3) read_fault_params(p.coord_file, off, csize, coords_buf);
        read_fault_params(p.vs_file, off, csize, vs_buf);
        read_fault_params(p.rho_file, off, csize, rho_buf);

        /* looping over each subfault */
        for (k=0; k<csize; k++) {
            calculate_moment_rate(time_buf, nt, p.source_dx, xx_buf[k], yy_buf[k], zz_buf[k], xz_buf[k], yz_buf[k], xy_buf[k],
                                    slip_buf[k], psv_buf[k], trup_buf[k], strike_buf[k], dip_buf[k], rake_buf[k], vs_buf[k], rho_buf[k], 
                                        rank, p.dc, p.median_ts, p.truptot, p.tp_psv_coef);

            
            /* filter if desired */
            if (p.filter) {
                if (rank==0) fprintf(stderr, "filtering source with parameters listed above.\n");
                xapiir_(xx_buf[k], &nt, p.aproto, &p.trbndw, &p.a, &p.iord, p.ftype, &p.hp, &p.lp, &p.dt, &p.npas);
                xapiir_(yy_buf[k], &nt, p.aproto, &p.trbndw, &p.a, &p.iord, p.ftype, &p.hp, &p.lp, &p.dt, &p.npas);
                xapiir_(zz_buf[k], &nt, p.aproto, &p.trbndw, &p.a, &p.iord, p.ftype, &p.hp, &p.lp, &p.dt, &p.npas);
                xapiir_(xz_buf[k], &nt, p.aproto, &p.trbndw, &p.a, &p.iord, p.ftype, &p.hp, &p.lp, &p.dt, &p.npas);
                xapiir_(yz_buf[k], &nt, p.aproto, &p.trbndw, &p.a, &p.iord, p.ftype, &p.hp, &p.lp, &p.dt, &p.npas);
                xapiir_(xy_buf[k], &nt, p.aproto, &p.trbndw, &p.a, &p.iord, p.ftype, &p.hp, &p.lp, &p.dt, &p.npas);
            }

            /* determine subfault location assuming fault is striking along x component */
            xil_buf[k] = p.x_start + ((s0 + k) % p.nx * source_sim_dx_rat); 

            // rough
            if (p.proj == 3) {
                yil_buf[k] = p.y_start + rint((coords_buf[k] - p.faultn_coord) / p.sim_dx);
            }
            // planar
            else {
                yil_buf[k] = p.y_start;
            }
            zil_buf[k] = p.z_start + ((s0 + k) / p.nx * source_sim_dx_rat);  

            /* compute local moment filtered source */
            moment += calculate_moment(nt, p.dt, xx_buf[k], yy_buf[k], zz_buf[k], xz_buf[k], yz_buf[k], xy_buf[k]);
            
        } // end loop subfault
            
        /* write moment rates */
        if (rank==master && debug) fprintf(stderr, "writing moment-rates to file.\n");
        write_momrate(p.momentrate_file, nt, p.nchunks, rank, csize, l, xil_buf, yil_buf, zil_buf, xx_buf, yy_buf, zz_buf, xz_buf, yz_buf, xy_buf);
        
    } // end main loop

    /* compute moment from all processes */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&moment, &global_moment, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
    if (rank == 0) {
        fprintf(stderr, "total moment: %e\n", global_moment);
        fprintf(stderr, "mw: %e\n", (2.0/3.0)*(log10(global_moment)-9.1));
    }

    /* free buffers */
    for (l=0; l<csize; l++) {
        free(xx_buf[l]);
        free(yy_buf[l]);
        free(zz_buf[l]);
        free(xz_buf[l]);
        free(yz_buf[l]);
        free(xy_buf[l]);
    }
    free(xx_buf);
    free(yy_buf);
    free(zz_buf);
    free(xz_buf);
    free(yz_buf);
    free(xy_buf);
    
    // 1d buffers
    free(psv_buf);
    free(trup_buf);
    free(strike_buf);
    free(dip_buf);
    free(rake_buf);
    free(coords_buf);
    free(slip_buf);
    free(vs_buf);
    free(rho_buf);
    free(xil_buf);
    free(yil_buf);
    free(zil_buf);
    free(time_buf);
    
    /* finalize mpi */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

/* FUNCTIONS */
void calculate_moment_rate(float *time, int nt, float dx, float *xx, float *yy, float *zz, float *xz, float *yz, float *xy, 
                                    float slip, float psv, float trup, float strike, float dip, float rake, 
                                        float vs, float rho, int rank, float dc, float median_ts, float truptot, float tp_psv_coef) {
    float *sliprate; 
    float momentrate;
    float ts, tr;
    float dt;
    int i;
    float srad, drad, rrad;
    // float ratio = 0.04235592747; // hard-coded assuming fs_max = 15.0 Hz
    //float dc = 0.11;
    //float median_ts = 1.0/15.0;
    float ratio;
    //float truptot=18.5;
    float temp_psv, treff;

    tr=0.0;
    ts=0.0;
    // compute dt
    dt = time[1] - time[0];
    // compute ratio of slip to vpeak
    ratio=slip/psv;
    // modify vpeak so ratio is bounded at 2.0
    temp_psv=psv;
    if (ratio > 2) {
        temp_psv = slip / 2;
    }
    // limit the vpeak to 0.1 m/s
    if (temp_psv < 0.1) {
        temp_psv = 0.1;
    }
    // compute ts using altered psv
    ts = tp_psv_coef * dc / temp_psv;
    // impose minimum value of 0.1 m/s to cap the ts at 0.825 seconds.
    if (ts >= 0.824999) {
        ts = median_ts;
    }
    treff = truptot-trup;
    // define 0 initialized vector to hold the slip-rate
    sliprate = zeros_f(nt);
    // for subfaults with abnormally large rupture times, assume they did not rupture
    if (treff >= 0.0) {
        tr=compute_tr(slip, treff);
        if (tr <= ts) {
            // arbitrary increment
            tr = ts + 100*dt;
        }
        // this should not happen, but its here just in case.
        else if (tr >= truptot) {
            tr = truptot; 
        }
        // compute tinti function using these parameters
        sliprate = tinti(time, nt, tr, ts, trup, slip);
    }
    // output debugging information
    fprintf(stderr, "(%d) psv: %f slip: %f trup: %f strike: %f dip: %f rake: %f ts: %f tr: %f nt: %i\n", rank, psv, slip, trup, strike, dip, rake, ts, tr, nt);


    /* DEPRECATED : Exponential Source-time function 
    tpeak = slip / (exp(1)*psv);

    for (i=0; i<nt; i++) {
        // calculate slip-rates
        if (time[i] > trup ) {
            sliprate = slip / tpeak * (time[i] - trup) / tpeak * exp(-((time[i] - trup) / tpeak));
        } 
        // prevent singularity when slip = 0; tpeak=0, and sliprate = inf
        else if (slip < 0.0001) {
            sliprate = 0.0;
        }
        else {
            sliprate = 0.0;
        }
        // muAD = moment
    */
        
    for (i=0; i<nt; i++) {
        // convert to moment-rate using Mo = uAD
        momentrate = sliprate[i] * vs*vs*rho * dx*dx;

        // convert to moment-rate tensor from strike dip and rake based on aki and richards using awp-odc sign conventions.
        srad = strike * M_PI / 180.0;
        drad = dip * M_PI / 180.0;
        rrad = rake * M_PI / 180.0;
        xx[i] = (momentrate * (sin(drad)*cos(rrad)*sin(2*srad) - sin(2*drad)*sin(rrad)*cos(srad)*cos(srad)));
        yy[i] = (-momentrate * (sin(drad)*cos(rrad)*sin(2*srad) + sin(2*drad)*sin(rrad)*sin(srad)*sin(srad)));
        zz[i] = (momentrate * (sin(2*drad)*sin(rrad)));
        xz[i] = (momentrate * (cos(drad)*cos(rrad)*sin(srad) - cos(2*drad)*sin(rrad)*cos(srad)));
        yz[i] = (momentrate * (cos(drad)*cos(rrad)*cos(srad)+cos(2*drad)*sin(rrad)*sin(srad)));
        xy[i] = (momentrate * (sin(drad)*cos(rrad)*cos(2*srad) + 0.5*sin(2*drad)*sin(rrad)*sin(2*srad)));
    }
    // release memory
    free(sliprate);
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

float compute_tr(float slip, float treff) {
    float tr;
    tr=2.596*slip+0.108*treff;
    return tr;
}
