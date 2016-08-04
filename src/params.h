typedef struct {
   int nchunks, nx, nz, nt;
   float dt, dx, rt;
   int x_start, y_start, z_start;
   float mean_faultn_coord; // remove fault start
   char psv_file[256], trup_file[256];
   char strike_file[256], dip_file[256], rake_file[256];
   char momentrate_file[256], moment_file[256], coord_file[256];
   int iord, npas;
   char aproto[6], ftype[6];
   float lp, hp, trbndw, a;
} param;

void read_settings(param *P, char *paramsfile);