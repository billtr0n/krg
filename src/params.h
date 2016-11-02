typedef struct {
   int nchunks, nx, nz, nt;
   float dt, source_dx, rt;
   float sim_dx;
   int x_start, y_start, z_start;
   float faultn_coord;
   char psv_file[256], trup_file[256];
   char strike_file[256], dip_file[256], rake_file[256];
   char rho_file[256], vs_file[256];
   char momentrate_file[256], slip_file[256], coord_file[256];
   int iord, npas;
   char aproto[6], ftype[6];
   float lp, hp, trbndw, a;
} param;

void read_settings(param *P, char *paramsfile);
