#include "io.h"
#include<vector>
#include<string>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<sys/types.h>
#include<sys/wait.h>
#include<unistd.h>
using namespace std;

static char encoding_table[64] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                                  'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
                                  'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
                                  'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
                                  'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
                                  'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
                                  'w', 'x', 'y', 'z', '0', '1', '2', '3',
                                  '4', '5', '6', '7', '8', '9', '+', '/'};
static int mod_table[64] = {0, 2, 1};

FILE * io[2];
int child_pid;

void error_message(const char * file, int line)
{
	char message[256];
	sprintf(message, "ERROR: %s:%d", file, line);
	perror(message);
	exit(-1);
}

void load_environ_file(int * nfname, char * filename, int * flen, 
		       char * fnam, int * clen, char * cnam, int * vlen, char * vnam)
{
	string fname;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}

	*clen = -1;
	*flen = -1;
	*vlen = -1;
	string line, header, section;
	ifstream in(fname.c_str());
	if ( ! in.good() ) {
		error_message(__FILE__, __LINE__);
	}
	getline(in, header);
	if ( header != "ENVIRONMENT FILE EXFEM" ) 
	{
		printf("ERROR: incorrect header.");
	}
	while( getline(in, section) )
	{
		if ( section == "TMP:" ) 
		{
			getline(in, line);
			*flen = std::min(int(line.size()), 256);
			for ( int i = 0; i < 256; i++ )
			{
				if ( i < *flen ) 
				{
					fnam[i] = line[i];
				} else
				{
					fnam[i] = ' ';
				}
			}
		}
		else if ( section == "VISUALIZATION:" ) 
		{
			getline(in, line);
			*vlen = std::min(int(line.size()), 256);
			for ( int i = 0; i < 256; i++ )
			{
				if ( i < *vlen ) 
				{
					vnam[i] = line[i];
				} else
				{
					vnam[i] = ' ';
				}
			}
		}
		else if ( section == "CVM PATH:" ) 
		{
			getline(in, line);
			*clen = std::min(int(line.size()), 256);
			for ( int i = 0; i < 256; i++ )
			{
				if ( i < *clen ) 
				{
					cnam[i] = line[i];
				} else
				{
					cnam[i] = ' ';
				}
			}
		}
		else if ( section != "" )
		{
			cout << "ERROR: Not recognized section: " << section << endl;
		}
	}
}

void load_problem_file(int * nfname, char * filename, double * height, 
		       double * width, double * time, double * inttime,
		       double * chper, double * nperlen, double * sperlen,
		       int * nmats, char ** cmats, double ** vmats,
		       char * prob, int * nppars, double ** vppars, 
		       char * geom, int * ngpars, double ** vgpars )
{
	string fname;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	string line, header, section;
	ifstream in(fname.c_str());
	if ( ! in.good() ) {
		error_message(__FILE__, __LINE__);
	}
	getline(in, header);
	if ( header != "PROBLEM FILE EXFEM" ) 
	{
		printf("ERROR: incorrect header.");
	}
	while( getline(in, section) ) {
		if ( section == "DOMAIN:" )
		{
			getline(in, line);
			istringstream is(line);
			is >> *width >> *height;
		}
		else if ( section == "CH PERIOD:" )
		{
			getline(in, line);
			istringstream is(line);
			is >> *chper;
		}
		else if ( section == "TIME:" )
		{
			getline(in, line);
			istringstream is(line);
			is >> *time >> *inttime;
		}
		else if ( section == "CONTROL:" ) 
		{
			getline(in, line);
			istringstream is(line);
			is >> *nperlen >> *sperlen;
		}
		else if ( section == "MATERIALS:" )
		{
			getline(in, line);
			istringstream is(line);
			is >> *nmats;
			*cmats = (char *)malloc((*nmats) * 2 * sizeof(char));
			*vmats = (double *)malloc((*nmats) * 3 * sizeof(double));
			for ( int i = 0; i < (*nmats); i++ )
			{
				string s;
				getline(in, line);
				istringstream iss(line);	
				iss >> s;
				(*cmats)[i*2+0] = s[0];
				iss >> (*vmats)[i*3+0];
				iss >> s;
				(*cmats)[i*2+1] = s[0];
				iss >> (*vmats)[i*3+1];
				iss >> (*vmats)[i*3+2];
			}
		}
		else if ( section == "PROBLEM:" )
		{
			getline(in, line);
			for ( int i = 0; i < 6; i++ )
			{
				if ( i < int(line.size()) )
				{
					prob[i] = line[i];
				}
				else
				{
					prob[i] = ' ';
				}
			}
			getline(in, line);
			istringstream is(line);
			is >> *nppars;
			*vppars = (double *)malloc((*nppars) * sizeof(double));
			for ( int i = 0; i < (*nppars); i++ )
			{
				is >> (*vppars)[i];
			}
		}
		else if ( section == "GEOMETRY:" )
		{
			getline(in, line);
			for ( int i = 0; i < 6; i++ )
			{
				if ( i < int(line.size()) )
				{
					geom[i] = line[i];
				}
				else
				{
					geom[i] = ' ';
				}
			}
			getline(in, line);
			istringstream is(line);
			is >> *ngpars;
			*vgpars = (double *)malloc((*ngpars) * sizeof(double));
			for ( int i = 0; i < (*ngpars); i++ )
			{
				is >> (*vgpars)[i];
			}
		}
	}
}

void save_rank_file(int * nfname, char * filename, int * myrank)
{
	string fname;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(), "rb+");
	int nread, nwritten;
	int rankpos;
	if ( file == NULL )
	{
		file = fopen(fname.c_str(), "wb+");
		if ( file == NULL ) 
			error_message(__FILE__, __LINE__);
		int zero = 0;
		nwritten = fwrite(&zero, sizeof(int), 1, file);
		if ( nwritten != 1 )
			error_message(__FILE__,__LINE__);
		nwritten = fwrite(&zero, sizeof(int), 1, file);
		if ( nwritten != 1 )
			error_message(__FILE__,__LINE__);
		fseek(file, 0, SEEK_SET);
	}
	nread = fread(&rankpos, sizeof(int), 1, file);
	if ( nread != 1 ) 
		error_message(__FILE__,__LINE__);
	fseek(file, (sizeof(int)*(2l+size_t(rankpos))), SEEK_SET);
	nwritten = fwrite(myrank, sizeof(int), 1, file);
	if ( nwritten != 1 ) 
		error_message(__FILE__,__LINE__);
	rankpos++;
	fseek(file, 0, SEEK_SET);
	nwritten = fwrite(&rankpos, sizeof(int), 1, file);
	if ( nwritten != 1 )
		error_message(__FILE__,__LINE__);
	fclose(file);
}

void restart_rank_file(int * nfname, char * filename)
{
	string fname;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	
	FILE * file = fopen(fname.c_str(), "rb+");
	int nwritten;
	int zero = 0;
	if ( file == NULL ) {
		cout << "Has a problem been solved in this machine?" << endl;
		error_message(__FILE__,__LINE__);
	}
	fseek(file, sizeof(int), SEEK_SET);
	nwritten = fwrite(&zero, sizeof(int), 1, file);
	if ( nwritten != 1 )
		error_message(__FILE__,__LINE__);
	fclose(file);
}

void load_rank_file(int * nfname, char * filename, int * nranks, int * myrank)
{
	string fname;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(), "rb+");
	if ( file == NULL )
		error_message(__FILE__,__LINE__);
	size_t nread, nwritten;
	int rankpos;
	nread = fread(nranks, sizeof(int), 1, file);
	if ( nread != 1 )
		error_message(__FILE__,__LINE__);
	nread = fread(&rankpos, sizeof(int), 1, file);
	if ( rankpos < *nranks ) {
		fseek(file, rankpos*sizeof(int), SEEK_CUR);
		nread = fread(myrank, sizeof(int), 1, file);
		if ( nread != 1 )
			error_message(__FILE__,__LINE__);
		fseek(file, sizeof(int), SEEK_SET);
		rankpos++;
		nwritten = fwrite(&rankpos, sizeof(int), 1, file);
		if ( nwritten != 1 )
			error_message(__FILE__,__LINE__);
	} else {
		cout << "ERROR: going beyond rankpos.\n More processors during visualization?" << endl;
	}
	fclose(file);
}

void remove_rank_file(int * nfname, char * filename)
{
	string fname;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	remove(fname.c_str());
}

void save_binmesh(int * nfname, char * filename, int * nelems, int * eptr,
	          int * eind, bool * ediv, int * nnodes, double * nodes)
{
	string fname;
	size_t nwritten;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(), "wb");
	if ( file == NULL )
	{
		cerr << "\"" << fname << "\"" << endl;
		error_message(__FILE__,__LINE__);
	}
	nwritten = fwrite(nelems, sizeof(int), size_t(1), file);
	if ( nwritten != 1 )
		error_message(__FILE__,__LINE__);
	nwritten = fwrite(eptr, sizeof(int), size_t( (*nelems) + 1), file);
	if ( nwritten != size_t((*nelems)+1) ) 
		error_message(__FILE__,__LINE__);
	nwritten = fwrite(eind, sizeof(int), size_t(eptr[*nelems]), file);
	if ( nwritten != size_t(eptr[*nelems]) )
		error_message(__FILE__,__LINE__);
	nwritten = fwrite(ediv, sizeof(bool), size_t(eptr[*nelems]), file);
	if ( nwritten != size_t(eptr[*nelems]) )
		error_message(__FILE__,__LINE__);
	nwritten = fwrite(nnodes, sizeof(int), size_t(1), file);
	if ( nwritten != 1 )
		error_message(__FILE__,__LINE__);
	nwritten = fwrite(nodes, sizeof(double), size_t( 2 * (*nnodes) ), file);
	if ( nwritten != size_t( 2 * (*nnodes) ) )
		error_message(__FILE__,__LINE__);
	fclose(file);
}

void load_binmesh(int * nfname, char * filename, 
	       int * nelems, int ** eptr, int ** eind, bool ** ediv,
	       int * nnodes, double ** nodes)
{
	string fname;
	size_t nread;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(), "rb");
	if ( file == NULL )
	{
		error_message(__FILE__,__LINE__);
	}
	nread = fread(nelems, sizeof(int), 1, file);
	if ( nread != size_t(1) )
		error_message(__FILE__,__LINE__);
	*eptr = (int *)malloc(sizeof(int)*size_t((*nelems)+1));
	nread = fread(*eptr, sizeof(int), size_t((*nelems)+1), file);
	if ( nread != size_t((*nelems)+1) )
		error_message(__FILE__,__LINE__);
	*eind = (int *)malloc(sizeof(int)*size_t((*eptr)[*nelems]));
	nread = fread(*eind, sizeof(int), size_t((*eptr)[*nelems]), file);
	if ( nread != size_t((*eptr)[*nelems]) )
		error_message(__FILE__,__LINE__);
	*ediv = (bool *)malloc(sizeof(bool)*size_t((*eptr)[*nelems]));
	nread = fread(*ediv, sizeof(bool), size_t((*eptr)[*nelems]), file);
	if ( nread != size_t((*eptr)[*nelems]) )
		error_message(__FILE__,__LINE__);
	nread = fread(nnodes, sizeof(int), size_t(1), file);
	if ( nread != 1 )
		error_message(__FILE__,__LINE__);
	*nodes = (double *)malloc(sizeof(double)*size_t( 2*(*nnodes) ));
	nread = fread(*nodes, sizeof(double), size_t( 2*(*nnodes) ), file);
	if ( nread != size_t( 2*(*nnodes) ) )
		error_message(__FILE__,__LINE__);
	fclose(file);
}

void save_binmass(int * nfname, char * filename, int * ndofs, double * mass)
{
	string fname;
	size_t nwritten;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(), "wb");
	if ( file == NULL )
	{
		error_message(__FILE__,__LINE__);
	}
	nwritten = fwrite(ndofs, sizeof(int), size_t(1), file);
	if ( nwritten != size_t(1) )
		error_message(__FILE__,__LINE__);
	nwritten = fwrite(mass, sizeof(double), size_t(*ndofs), file);
	if ( nwritten != size_t(*ndofs) )
		error_message(__FILE__,__LINE__);
	fclose(file);
}

void load_binmass(int * nfname, char * filename, int * ndofs, double ** mass)
{
	string fname;
	size_t nread;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(), "rb");
	if ( file == NULL )
	{
		error_message(__FILE__,__LINE__);
	}
	nread = fread(ndofs, sizeof(int), size_t(1), file);
	if ( nread != size_t(1) )
		error_message(__FILE__,__LINE__);
	*mass = (double *)malloc(sizeof(double)*size_t(*ndofs));
	nread = fread(*mass, sizeof(double), size_t(*ndofs), file);
	if ( nread != size_t(*ndofs) )
		error_message(__FILE__,__LINE__);
	fclose(file);
}

void save_binstep(int * nfname, char * filename, char * mode, 
		   int * step, int * ndofs, double * u)
{
	string fname;
	size_t nwritten;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file;
        if ( mode[0] == 'n' ) 
		file = fopen(fname.c_str(),"wb");
	else
		file = fopen(fname.c_str(),"ab");
	if ( file == NULL )
	{
		error_message(__FILE__,__LINE__);
	}
	nwritten = fwrite(step, sizeof(int), 1, file);
	if ( nwritten != 1 )
		error_message(__FILE__,__LINE__);
	nwritten = fwrite(ndofs, sizeof(int), 1, file);
	if ( nwritten != 1 ) 
		error_message(__FILE__,__LINE__);
	nwritten = fwrite(u, sizeof(double), size_t(*ndofs), file);
	if ( nwritten != size_t(*ndofs) ) 
		error_message(__FILE__,__LINE__);
	fclose(file);
}

void load_binstep(int * nfname, char * filename, int * step, int * ndofs, double ** u, int * status)
{
	string fname;
	size_t nread;
	int fstep;
	int rstep;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(),"rb");
	if ( file == NULL )
	{
		error_message(__FILE__,__LINE__);
	}
	nread = fread(&fstep, sizeof(int), 1, file);
	if ( nread != 1 )
		error_message(__FILE__,__LINE__);
	rstep = (*step)+fstep;
	nread = fread(ndofs, sizeof(int), 1, file);
	if ( nread != 1 ) 
		error_message(__FILE__,__LINE__);
	fseek(file, 0l, SEEK_END);
	size_t end = ftell(file);
	size_t jmp = size_t(*step)*((sizeof(int)*2l)+(size_t(*ndofs)*sizeof(double)));
	if ( jmp >= end ) 
	{
		*status = -1;
		return;
	}
	fseek(file, jmp, SEEK_SET);
	nread = fread(&fstep, sizeof(int), 1, file);
	if ( nread != 1 )
		error_message(__FILE__,__LINE__);
	if ( fstep != rstep )
		cout << "Error: reading wrong step" << endl;
	nread = fread(ndofs, sizeof(int), 1, file);
	if ( nread != 1 ) 
		error_message(__FILE__,__LINE__);
	*u = (double *)malloc(sizeof(double)*size_t(*ndofs));
	nread = fread(*u, sizeof(double), size_t(*ndofs), file);
	if ( nread != size_t(*ndofs) ) 
		error_message(__FILE__,__LINE__);
	fclose(file);
	*status = 0;
}

void visual_params(char * type, double * params)
{
	string stype;
	cin >> stype;
	for ( size_t i = 0; i < 6; i++ )
	{
		if ( i < stype.size() )
		{
			type[i] = stype[i];
		}
		else
		{
			type[i] = ' ';
		}
	}
	if ( strncmp(type, "SHEET", 5) == 0 ) {
		cout << "Parameters for sheet: xmin, xmax, dist\n";
		cin >> params[0] >> params[1] >> params[2];
	} else if ( strncmp(type, "ANIMAT", 6) == 0 ) {
		cout << "Parameters for animation: xmin, xmax, ymin, ymax, dist\n";
		cin >> params[0] >> params[1] >> params[2] >> params[3] >> params[4];
	} else {
		cerr << "ERROR: The type of output is not recognized.\n";
	}
}

void open_cvm_pipe(int * npath, char * path)
{
	int in[2];
	if ( pipe(in) == -1 ) 
		error_message(__FILE__,__LINE__);
	int out[2];
	if ( pipe(out) == -1 )
		error_message(__FILE__,__LINE__);
	child_pid = fork();
	if ( child_pid == 0 )
	{
		char cpath[256];
		int i;
		for ( i = 0; (i < *npath && i < 256); i++ )
		{
			cpath[i] = path[i];
		}
		cpath[i] = '\0';
		if ( close(in[1]) != 0 ) 
			error_message(__FILE__,__LINE__);
		if ( close(out[0]) != 0 ) 
			error_message(__FILE__,__LINE__);
		if ( dup2(in[0], fileno(stdin)) == -1 )
			error_message(__FILE__,__LINE__);
		if ( dup2(out[1], fileno(stdout)) == -1 )
			error_message(__FILE__,__LINE__);
		if ( chdir(cpath) != 0 ) 
			error_message(__FILE__,__LINE__);
		if ( execlp("./slicerecon", "./slicerecon", "cvm", (char *)NULL) == -1 )
			error_message(__FILE__,__LINE__);
	}
	else
	{
		if ( close(in[0]) != 0 ) 
			error_message(__FILE__,__LINE__);
		if ( close(out[1]) != 0 )
			error_message(__FILE__,__LINE__);
		io[0] = fdopen(in[1], "w");
		if ( io[0] == NULL )
			error_message(__FILE__,__LINE__);
		io[1] = fdopen(out[0], "r");
		if ( io[1] == NULL )
			error_message(__FILE__,__LINE__);	
	}
}

void cvm_mat(double * point, int * mat )
{
	fprintf(io[0], "%lf %lf %lf\n", point[0], point[1], point[2]);
	fflush(io[0]);
	int i = fscanf(io[1], "%d", mat);
	if ( i != 1 ) 
	{
		error_message(__FILE__,__LINE__);
	}
}

void close_cvm_pipe()
{
	int status;
	fprintf(io[0], "end\n");
	fflush(io[0]);
	child_pid = waitpid(child_pid, &status, 0);
	fclose(io[0]);
	fclose(io[1]);
}

void base64_encode(int nin, const unsigned char *in, int & nout, char ** out) {
	nout = 4 * ((nin + 2) / 3);
	(*out) = (char *)malloc(sizeof(char)*(nout+1));
	if ((*out) == NULL)
	{
		error_message(__FILE__,__LINE__);
	}
	for (int i = 0, j = 0; i < nin;) {
		unsigned int octet_a = i < nin ? in[i++] : 0;
		unsigned int octet_b = i < nin ? in[i++] : 0;
		unsigned int octet_c = i < nin ? in[i++] : 0;
		unsigned int triple = (octet_a << 0x10) + (octet_b << 0x08) + octet_c;
		(*out)[j++] = encoding_table[(triple >> 3 * 6) & 0x3F];
		(*out)[j++] = encoding_table[(triple >> 2 * 6) & 0x3F];
		(*out)[j++] = encoding_table[(triple >> 1 * 6) & 0x3F];
		(*out)[j++] = encoding_table[(triple >> 0 * 6) & 0x3F];
	}

	for (int i = 0; i < mod_table[nin % 3]; i++)
		(*out)[nout - 1 - i] = '=';
	(*out)[nout]='\0';
}

void num_to_base64(int num, char out[9])
{
	unsigned char * ch = (unsigned char *)&num;
	unsigned int octet_a = ch[0];
	unsigned int octet_b = ch[1];
	unsigned int octet_c = ch[2];
	unsigned int triple = (octet_a << 0x10) + (octet_b << 0x08) + octet_c;
	out[0] = encoding_table[(triple >> 3 * 6) & 0x3F];
	out[1] = encoding_table[(triple >> 2 * 6) & 0x3F];
	out[2] = encoding_table[(triple >> 1 * 6) & 0x3F];
	out[3] = encoding_table[(triple >> 0 * 6) & 0x3F];
	octet_a = ch[4];
	triple = (octet_a << 0x10);
	out[4] = encoding_table[(triple >> 3 * 6) & 0x3F];
	out[5] = encoding_table[(triple >> 2 * 6) & 0x3F];
	out[6] = '=';
	out[7] = '=';
	out[8] = '\0';
}

void mesh_to_base64(int nelems, int * eptr, int * eind, int& ncout, char ** cout, int& noout, char ** oout )
{
	int nconn = eptr[nelems]-1;
	int * imesh = (int *)malloc(nconn*sizeof(int));
	int * ioff  = (int *)malloc(nelems*sizeof(int));
	for ( int i = 0; i < nelems; i++ )
	{
		ioff[i] = (eptr[i+1]-1);
		for ( int j = eptr[i]-1; j < eptr[i+1]-1; j++ )
		{
			imesh[j] = (eind[j]-1);
		}
	}
	base64_encode(nelems*sizeof(int), (unsigned char *)ioff, noout, oout);
	base64_encode(nconn*sizeof(int), (unsigned char *)imesh, ncout, cout);
	free(ioff);
	free(imesh);
}

void pointdata_to_base64(int npoints, double * points, int& nout, char ** out)
{
	float * fpoints = (float *)malloc(sizeof(float)*3*npoints);
	for ( int i = 0; i < npoints; i++ )
	{
		fpoints[3*i]   = points[2*i];
		fpoints[3*i+1] = points[2*i+1];
		fpoints[3*i+2] = 0.0;
	}
	base64_encode(3*sizeof(float)*npoints, (unsigned char *)fpoints,  nout, out);
	free(fpoints);
}

void write_vtk_unst_field(int * nfname, char * filename, int * nelems, 
			  int * ptr, int * ind, int * npoints,
			  double * points, int * nfieldname, 
			  char * fieldname, double * field)
{
	string fname, fldname;
	int nwritten, n1, n2;
	int offsets[6];
	char * b64text1, * b64text2;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	for ( int i = 0; i < *nfieldname; i++ )
	{
		fldname += fieldname[i];
	}
	FILE * file = fopen(fname.c_str(),"w");
	if ( file == NULL ) 
	{
		error_message(__FILE__,__LINE__);
	}
	offsets[0] = 0;
	offsets[1] = offsets[0] + 4 * ((((*npoints)*sizeof(float)*3) + 2) / 3) + 8;
	offsets[2] = offsets[1] + 4 * ((((*npoints)*sizeof(float)*3) + 2) / 3) + 8;
	offsets[3] = offsets[2] + 4 * ((((ptr[*nelems]-1)*sizeof(int)) + 2) / 3) + 8;
	offsets[4] = offsets[3] + 4 * ((((*nelems)*sizeof(int)) + 2) / 3) + 8;
	offsets[5] = offsets[4] + 4 * ((((*nelems)*sizeof(int)) + 2) / 3) + 8;
	
        nwritten = fprintf(file, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "<UnstructuredGrid>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", *npoints, *nelems);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	
	nwritten = fprintf(file, "<PointData>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "<DataArray type=\"Float32\" Name=\"%s\"", fldname.c_str());
	nwritten = fprintf(file, " NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\">\n", offsets[0]);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "\n</DataArray>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "</PointData>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	
	nwritten = fprintf(file, "<CellData>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "</CellData>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	
	nwritten = fprintf(file, "<Points>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "<DataArray type=\"Float32\" Name=\"points\"");
	nwritten = fprintf(file, " NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\">\n", offsets[1]);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "\n</DataArray>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "</Points>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);

	nwritten = fprintf(file, "<Cells>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	
	nwritten = fprintf(file, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%d\">\n", offsets[2]);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "\n</DataArray>\n");
	
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "<DataArray type=\"Int32\" Name=\"offsets\" Format=\"appended\" offset=\"%d\">\n", offsets[3]);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "\n</DataArray>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	
	nwritten = fprintf(file, "<DataArray type=\"Int32\" Name=\"types\" Format=\"appended\" offset=\"%d\">\n", offsets[4]);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "\n</DataArray>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);

	nwritten = fprintf(file, "</Cells>\n");
	
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "</Piece>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "</UnstructuredGrid>\n\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	
	nwritten = fprintf(file, "<AppendedData encoding=\"base64\">\n\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	
	fputc('_', file);
	
	char out[9];
	
	
	// printing field.
	pointdata_to_base64(*npoints, field, n1, &b64text1);
	num_to_base64(n1, out);
	nwritten = fputs(out, file);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fputs(b64text1, file);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	free(b64text1);
	
	// printing points.
	pointdata_to_base64(*npoints, points, n1, &b64text1);
	num_to_base64(n1, out);
	nwritten = fputs(out, file);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fputs(b64text1, file);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	free(b64text1);

	// printing mesh.
	mesh_to_base64(*nelems, ptr, ind, n1, &b64text1, n2, &b64text2);
	num_to_base64(n1, out);
	nwritten = fputs(out, file);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fputs(b64text1, file);
	free(b64text1);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	
	num_to_base64(n2, out);
	nwritten = fputs(out, file);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fputs(b64text2, file);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	free(b64text2);
	
	nwritten = fputs(out, file);
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	for ( int i = 0; i < *nelems; i+=3 )
		nwritten = fputs("BwAAAAcAAAAHAAAA", file);
	if ( (*nelems)%3 == 1 )
		nwritten = fputs("BwAAAA==", file);
	else if ( (*nelems)%3 == 2 )
		nwritten = fputs("BwAAAAcAAAA=", file);

	nwritten = fprintf(file, "\n\n</AppendedData>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	nwritten = fprintf(file, "</VTKFile>\n");
	if ( nwritten < 0 ) error_message(__FILE__,__LINE__);
	fclose(file);
}

void write_sheet_points(int * nfname, char * filename, int * n, double * points)
{
	string fname;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(),"w");
	fprintf(file, "POINTS = [ ");
	for ( int i = 0; i < *n; i++ ) 
	{
		fprintf(file, "%.10lf ", points[i]);
	}
	fprintf(file, " ];");
	fprintf(file, "\nFIELD = [ ");
	fclose(file);
}

void write_sheet_step(int * nfname, char * filename, int * n, double * field)
{
	string fname;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(),"a");
	for ( int i = 0; i < *n; i++ ) 
	{
		fprintf(file, "%lf ", field[i]);
	}
	fprintf(file, "\n");
	fclose(file);
}

void write_sheet_close(int * nfname, char * filename, double * step)
{
	string fname;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(),"a");
	fprintf(file, " ];\n");
	fprintf(file, "SAMP = %lf ;", *step);
	fclose(file);
}

void save_bintparams(int * nfname, char * filename, double * deltat, int * iointerval)
{
	string fname;
	size_t nwritten;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(), "wb");
	if ( file == NULL )
	{
		error_message(__FILE__,__LINE__);
	}
	nwritten = fwrite(deltat, sizeof(double), size_t(1), file);
	if ( nwritten != size_t(1) )
		error_message(__FILE__,__LINE__);
	nwritten = fwrite(iointerval, sizeof(int), size_t(1), file);
	if ( nwritten != size_t(1) )
		error_message(__FILE__,__LINE__);
	fclose(file);
}

void load_bintparams(int * nfname, char * filename, double * deltat, int * iointerval)
{
	string fname;
	size_t nread;
	for ( int i = 0; i < *nfname; i++ )
	{
		fname += filename[i];
	}
	FILE * file = fopen(fname.c_str(), "rb");
	if ( file == NULL )
	{
		error_message(__FILE__,__LINE__);
	}
	nread = fread(deltat, sizeof(double), size_t(1), file);
	if ( nread != size_t(1) )
		error_message(__FILE__,__LINE__);
	nread = fread(iointerval, sizeof(int), size_t(1), file);
	if ( nread != size_t(1) )
		error_message(__FILE__,__LINE__);
	fclose(file);
}

