extern "C" {
#include<stdbool.h>

void load_environ_file(int * nfname, char * filename, int * flen, char * fnam, 
		       int * clen, char * cnam, int * vlen, char * vnam);
void load_problem_file(int * nfname, char * filename, double * height,
		       double * width, double * time, double * inttime,
		       double * chper, double * nperlen, double * sperlen,
		       int * nmats, char ** cmats, double ** vmats,
		       char * prob, int * nppars, double ** vppars,
		       char * geom, int * ngpars, double ** vgpars);
void save_binmesh(int * nfname, char * filename, int * nelems, int * eptr, int * eind, bool * ediv,
	          int * nnodes, double * nodes);
void load_binmesh(int * nfname, char * filename,  int * nelems, int ** eptr, int ** eind, bool ** ediv,
	          int * nnodes, double ** nodes);
void save_binmass(int * nfname, char * filename, int * ndofs, double * mass);
void load_binmass(int * nfname, char * filename, int * ndofs, double ** mass);
void save_binstep(int * nfname, char * filename, char * mode, int * step, int * ndofs, double * u);
void load_binstep(int * nfname, char * filename, int * step, int * ndofs, double ** u, int * status);

void visual_params(char * type, double * params);  
void remove_rank_file(int * nfname, char * filename);
void load_rank_file(int * nfname, char * filename, int * nranks, int * myrank);
void restart_rank_file(int * nfname, char * filename);
void save_rank_file(int * nfname, char * filename, int * myrank);
void open_cvm_pipe(int * npath, char * path);
void cvm_mat(double * point, int * mat);
void close_cvm_pipe();
void write_vtk_unst_field(int * nfname, char * filename, int * nelems, 
			  int * ptr, int * ind, int * npoints,
			  double * points, int * nfieldname, 
			  char * fieldname, double * field);
void write_sheet_close(int * nfn, char * fname, double * step);
void write_sheet_step(int * nfn, char * fname, int * n, double * field);
void write_sheet_points(int * nfn, char * fname, int * n, double * points);
void save_bintparams(int * nfname, char * filename, double * deltat, int * iointerval);
void load_bintparams(int * nfname, char * filename, double * deltat, int * iointerval);
};

