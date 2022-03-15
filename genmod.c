/* This function was automatically generated by CasADi */
#ifdef __cplusplus
extern "C" {
#endif

#ifdef CODEGEN_PREFIX
  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)
  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else /* CODEGEN_PREFIX */
  #define CASADI_PREFIX(ID) genmod_ ## ID
#endif /* CODEGEN_PREFIX */

#include <string.h>
#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif
#include <math.h>

#ifndef real_t
#define real_t double
#endif /* real_t */

#define to_double(x) (double) x
#define to_int(x) (int) x
/* Pre-c99 compatibility */
#if __STDC_VERSION__ < 199901L
real_t CASADI_PREFIX(fmin)(real_t x, real_t y) { return x<y ? x : y;}
#define fmin(x,y) CASADI_PREFIX(fmin)(x,y)
real_t CASADI_PREFIX(fmax)(real_t x, real_t y) { return x>y ? x : y;}
#define fmax(x,y) CASADI_PREFIX(fmax)(x,y)
#endif

#ifdef MATLAB_MEX_FILE
#define PRINTF mexPrintf
#else
#define PRINTF printf
#endif
real_t CASADI_PREFIX(sq)(real_t x) { return x*x;}
#define sq(x) CASADI_PREFIX(sq)(x)

real_t CASADI_PREFIX(sign)(real_t x) { return x<0 ? -1 : x>0 ? 1 : x;}
#define sign(x) CASADI_PREFIX(sign)(x)

void CASADI_PREFIX(fill)(real_t* x, int n, real_t alpha) {
  int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = alpha;
  }
}
#define fill(x, n, alpha) CASADI_PREFIX(fill)(x, n, alpha)


#ifdef MATLAB_MEX_FILE
real_t* CASADI_PREFIX(from_mex)(const mxArray *p, real_t* y, const int* sp, real_t* w) {
  if (!mxIsDouble(p) || mxGetNumberOfDimensions(p)!=2)
    mexErrMsgIdAndTxt("Casadi:RuntimeError","\"from_mex\" failed: Not a two-dimensional matrix of double precision.");
  int nrow = *sp++, ncol = *sp++, nnz = sp[ncol];
  const int *colind=sp, *row=sp+ncol+1;
  size_t p_nrow = mxGetM(p), p_ncol = mxGetN(p);
  const double* p_data = (const double*)mxGetData(p);
  bool is_sparse = mxIsSparse(p);
  mwIndex *Jc = is_sparse ? mxGetJc(p) : 0;
  mwIndex *Ir = is_sparse ? mxGetIr(p) : 0;
  if (p_nrow==1 && p_ncol==1) {
    double v = is_sparse && Jc[1]==0 ? 0 : *p_data;
    fill(y, nnz, v);
  } else {
    bool tr = false;
    if (nrow!=p_nrow || ncol!=p_ncol) {
      tr = nrow==p_ncol && ncol==p_nrow && (nrow==1 || ncol==1);
      if (!tr) mexErrMsgIdAndTxt("Casadi:RuntimeError","\"from_mex\" failed: Dimension mismatch.");
    }
    int r,c,k;
    if (is_sparse) {
      if (tr) {
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]+c*nrow]=0;
        for (c=0; c<p_ncol; ++c)
          for (k=Jc[c]; k<Jc[c+1]; ++k) w[c+Ir[k]*p_ncol] = p_data[k];
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) y[k] = w[row[k]+c*nrow];
      } else {
        for (c=0; c<ncol; ++c) {
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]]=0;
          for (k=Jc[c]; k<Jc[c+1]; ++k) w[Ir[k]]=p_data[k];
          for (k=colind[c]; k<colind[c+1]; ++k) y[k]=w[row[k]];
        }
      }
    } else {
      for (c=0; c<ncol; ++c) {
        for (k=colind[c]; k<colind[c+1]; ++k) {
          y[k] = p_data[row[k]+c*nrow];
        }
      }
    }
  }
  return y;
}
#define from_mex(p, y, sp, w) CASADI_PREFIX(from_mex)(p, y, sp, w)
#endif

#ifdef MATLAB_MEX_FILE
mxArray* CASADI_PREFIX(to_mex)(const int* sp, const real_t* x) {
  int nrow = *sp++, ncol = *sp++, nnz = sp[ncol];
  mxArray* p = mxCreateSparse(nrow, ncol, nnz, mxREAL);
  int i;
  mwIndex* j;
  for (i=0, j=mxGetJc(p); i<=ncol; ++i) *j++ = *sp++;
  for (i=0, j=mxGetIr(p); i<nnz; ++i) *j++ = *sp++;
  if (x) {
    double* d = (double*)mxGetData(p);
    for (i=0; i<nnz; ++i) *d++ = to_double(*x++);
  }
  return p;
}
#define to_mex(sp, x) CASADI_PREFIX(to_mex)(sp, x)
#endif

static const int CASADI_PREFIX(s0)[5] = {1, 1, 0, 1, 0};
#define s0 CASADI_PREFIX(s0)
static const int CASADI_PREFIX(s1)[6] = {2, 1, 0, 2, 0, 1};
#define s1 CASADI_PREFIX(s1)
static const int CASADI_PREFIX(s2)[8] = {2, 2, 0, 1, 3, 0, 0, 1};
#define s2 CASADI_PREFIX(s2)
static const int CASADI_PREFIX(s3)[7] = {2, 2, 0, 1, 2, 0, 1};
#define s3 CASADI_PREFIX(s3)
static const int CASADI_PREFIX(s4)[5] = {2, 1, 0, 1, 1};
#define s4 CASADI_PREFIX(s4)
static const int CASADI_PREFIX(s5)[5] = {4, 1, 0, 1, 3};
#define s5 CASADI_PREFIX(s5)
static const int CASADI_PREFIX(s6)[6] = {1, 2, 0, 1, 1, 0};
#define s6 CASADI_PREFIX(s6)
static const int CASADI_PREFIX(s7)[4] = {1, 1, 0, 0};
#define s7 CASADI_PREFIX(s7)
static const int CASADI_PREFIX(s8)[4] = {2, 1, 0, 0};
#define s8 CASADI_PREFIX(s8)
/* fd */
int fd(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=9.0000000000000002e-001;
  real_t a1=arg[1] ? arg[1][0] : 0;
  a1=(a0*a1);
  real_t a2=arg[1] ? arg[1][1] : 0;
  a1=(a1*a2);
  real_t a3=arg[4] ? arg[4][0] : 0;
  a1=(a1+a3);
  if (res[0]!=0) res[0][0]=a1;
  a1=1.0000000000000001e-001;
  a3=arg[3] ? arg[3][0] : 0;
  a1=(a1*a3);
  a0=(a0+a1);
  a0=(a0*a2);
  a2=arg[4] ? arg[4][1] : 0;
  a0=(a0+a2);
  if (res[0]!=0) res[0][1]=a0;
  return 0;
}

void fd_incref(void) {
}

void fd_decref(void) {
}

int fd_n_in(void) { return 5;}

int fd_n_out(void) { return 1;}

const char* fd_name_in(int i){
  switch (i) {
  case 0: return "t";
  case 1: return "x";
  case 2: return "u";
  case 3: return "th";
  case 4: return "w";
  default: return 0;
  }
}

const char* fd_name_out(int i){
  switch (i) {
  case 0: return "fd";
  default: return 0;
  }
}

const int* fd_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  case 2: return s0;
  case 3: return s0;
  case 4: return s1;
  default: return 0;
  }
}

const int* fd_sparsity_out(int i) {
  switch (i) {
  case 0: return s1;
  default: return 0;
  }
}

int fd_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 4;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_fd(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  int i, j;
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"fd\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"fd\" failed. Too many output arguments (%d, max 1)", resc);
  int *iw = 0;
  real_t w[13];
  const real_t* arg[5] = {0};
  if (--argc>=0) arg[0] = from_mex(argv[0], w, s0, w+9);
  if (--argc>=0) arg[1] = from_mex(argv[1], w+1, s1, w+9);
  if (--argc>=0) arg[2] = from_mex(argv[2], w+3, s0, w+9);
  if (--argc>=0) arg[3] = from_mex(argv[3], w+4, s0, w+9);
  if (--argc>=0) arg[4] = from_mex(argv[4], w+5, s1, w+9);
  real_t* res[1] = {0};
  --resc;
  res[0] = w+7;
  i = fd(arg, res, iw, w+9, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"fd\" failed.");
  if (res[0]) resv[0] = to_mex(s1, res[0]);
}
#endif

/* dfddx */
int dfddx(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=9.0000000000000002e-001;
  real_t a1=arg[1] ? arg[1][1] : 0;
  a1=(a0*a1);
  if (res[0]!=0) res[0][0]=a1;
  a1=arg[1] ? arg[1][0] : 0;
  a1=(a0*a1);
  if (res[0]!=0) res[0][1]=a1;
  a1=1.0000000000000001e-001;
  real_t a2=arg[3] ? arg[3][0] : 0;
  a1=(a1*a2);
  a0=(a0+a1);
  if (res[0]!=0) res[0][2]=a0;
  return 0;
}

void dfddx_incref(void) {
}

void dfddx_decref(void) {
}

int dfddx_n_in(void) { return 5;}

int dfddx_n_out(void) { return 1;}

const char* dfddx_name_in(int i){
  switch (i) {
  case 0: return "t";
  case 1: return "x";
  case 2: return "u";
  case 3: return "th";
  case 4: return "w";
  default: return 0;
  }
}

const char* dfddx_name_out(int i){
  switch (i) {
  case 0: return "dfddx";
  default: return 0;
  }
}

const int* dfddx_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  case 2: return s0;
  case 3: return s0;
  case 4: return s1;
  default: return 0;
  }
}

const int* dfddx_sparsity_out(int i) {
  switch (i) {
  case 0: return s2;
  default: return 0;
  }
}

int dfddx_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 3;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_dfddx(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  int i, j;
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dfddx\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dfddx\" failed. Too many output arguments (%d, max 1)", resc);
  int *iw = 0;
  real_t w[13];
  const real_t* arg[5] = {0};
  if (--argc>=0) arg[0] = from_mex(argv[0], w, s0, w+10);
  if (--argc>=0) arg[1] = from_mex(argv[1], w+1, s1, w+10);
  if (--argc>=0) arg[2] = from_mex(argv[2], w+3, s0, w+10);
  if (--argc>=0) arg[3] = from_mex(argv[3], w+4, s0, w+10);
  if (--argc>=0) arg[4] = from_mex(argv[4], w+5, s1, w+10);
  real_t* res[1] = {0};
  --resc;
  res[0] = w+7;
  i = dfddx(arg, res, iw, w+10, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dfddx\" failed.");
  if (res[0]) resv[0] = to_mex(s2, res[0]);
}
#endif

/* dfddw */
int dfddw(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=1.;
  if (res[0]!=0) res[0][0]=a0;
  if (res[0]!=0) res[0][1]=a0;
  return 0;
}

void dfddw_incref(void) {
}

void dfddw_decref(void) {
}

int dfddw_n_in(void) { return 5;}

int dfddw_n_out(void) { return 1;}

const char* dfddw_name_in(int i){
  switch (i) {
  case 0: return "t";
  case 1: return "x";
  case 2: return "u";
  case 3: return "th";
  case 4: return "w";
  default: return 0;
  }
}

const char* dfddw_name_out(int i){
  switch (i) {
  case 0: return "dfddw";
  default: return 0;
  }
}

const int* dfddw_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  case 2: return s0;
  case 3: return s0;
  case 4: return s1;
  default: return 0;
  }
}

const int* dfddw_sparsity_out(int i) {
  switch (i) {
  case 0: return s3;
  default: return 0;
  }
}

int dfddw_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 1;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_dfddw(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  int i, j;
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dfddw\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dfddw\" failed. Too many output arguments (%d, max 1)", resc);
  int *iw = 0;
  real_t w[11];
  const real_t* arg[5] = {0};
  if (--argc>=0) arg[0] = from_mex(argv[0], w, s0, w+9);
  if (--argc>=0) arg[1] = from_mex(argv[1], w+1, s1, w+9);
  if (--argc>=0) arg[2] = from_mex(argv[2], w+3, s0, w+9);
  if (--argc>=0) arg[3] = from_mex(argv[3], w+4, s0, w+9);
  if (--argc>=0) arg[4] = from_mex(argv[4], w+5, s1, w+9);
  real_t* res[1] = {0};
  --resc;
  res[0] = w+7;
  i = dfddw(arg, res, iw, w+9, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dfddw\" failed.");
  if (res[0]) resv[0] = to_mex(s3, res[0]);
}
#endif

/* dfddth */
int dfddth(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=1.0000000000000001e-001;
  real_t a1=arg[1] ? arg[1][1] : 0;
  a0=(a0*a1);
  if (res[0]!=0) res[0][0]=a0;
  return 0;
}

void dfddth_incref(void) {
}

void dfddth_decref(void) {
}

int dfddth_n_in(void) { return 5;}

int dfddth_n_out(void) { return 1;}

const char* dfddth_name_in(int i){
  switch (i) {
  case 0: return "t";
  case 1: return "x";
  case 2: return "u";
  case 3: return "th";
  case 4: return "w";
  default: return 0;
  }
}

const char* dfddth_name_out(int i){
  switch (i) {
  case 0: return "dfddth";
  default: return 0;
  }
}

const int* dfddth_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  case 2: return s0;
  case 3: return s0;
  case 4: return s1;
  default: return 0;
  }
}

const int* dfddth_sparsity_out(int i) {
  switch (i) {
  case 0: return s4;
  default: return 0;
  }
}

int dfddth_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 2;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_dfddth(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  int i, j;
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dfddth\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dfddth\" failed. Too many output arguments (%d, max 1)", resc);
  int *iw = 0;
  real_t w[10];
  const real_t* arg[5] = {0};
  if (--argc>=0) arg[0] = from_mex(argv[0], w, s0, w+8);
  if (--argc>=0) arg[1] = from_mex(argv[1], w+1, s1, w+8);
  if (--argc>=0) arg[2] = from_mex(argv[2], w+3, s0, w+8);
  if (--argc>=0) arg[3] = from_mex(argv[3], w+4, s0, w+8);
  if (--argc>=0) arg[4] = from_mex(argv[4], w+5, s1, w+8);
  real_t* res[1] = {0};
  --resc;
  res[0] = w+7;
  i = dfddth(arg, res, iw, w+8, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dfddth\" failed.");
  if (res[0]) resv[0] = to_mex(s4, res[0]);
}
#endif

/* ddfddxdth */
int ddfddxdth(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=1.0000000000000001e-001;
  if (res[0]!=0) res[0][0]=a0;
  return 0;
}

void ddfddxdth_incref(void) {
}

void ddfddxdth_decref(void) {
}

int ddfddxdth_n_in(void) { return 5;}

int ddfddxdth_n_out(void) { return 1;}

const char* ddfddxdth_name_in(int i){
  switch (i) {
  case 0: return "t";
  case 1: return "x";
  case 2: return "u";
  case 3: return "th";
  case 4: return "w";
  default: return 0;
  }
}

const char* ddfddxdth_name_out(int i){
  switch (i) {
  case 0: return "ddfddxdth";
  default: return 0;
  }
}

const int* ddfddxdth_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  case 2: return s0;
  case 3: return s0;
  case 4: return s1;
  default: return 0;
  }
}

const int* ddfddxdth_sparsity_out(int i) {
  switch (i) {
  case 0: return s5;
  default: return 0;
  }
}

int ddfddxdth_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 1;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_ddfddxdth(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  int i, j;
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"ddfddxdth\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"ddfddxdth\" failed. Too many output arguments (%d, max 1)", resc);
  int *iw = 0;
  real_t w[10];
  const real_t* arg[5] = {0};
  if (--argc>=0) arg[0] = from_mex(argv[0], w, s0, w+8);
  if (--argc>=0) arg[1] = from_mex(argv[1], w+1, s1, w+8);
  if (--argc>=0) arg[2] = from_mex(argv[2], w+3, s0, w+8);
  if (--argc>=0) arg[3] = from_mex(argv[3], w+4, s0, w+8);
  if (--argc>=0) arg[4] = from_mex(argv[4], w+5, s1, w+8);
  real_t* res[1] = {0};
  --resc;
  res[0] = w+7;
  i = ddfddxdth(arg, res, iw, w+8, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"ddfddxdth\" failed.");
  if (res[0]) resv[0] = to_mex(s5, res[0]);
}
#endif

/* g */
int g(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=arg[1] ? arg[1][0] : 0;
  a0=sq(a0);
  real_t a1=2.0000000000000001e-001;
  a1=(a1*a0);
  a0=arg[4] ? arg[4][0] : 0;
  a1=(a1+a0);
  if (res[0]!=0) res[0][0]=a1;
  return 0;
}

void g_incref(void) {
}

void g_decref(void) {
}

int g_n_in(void) { return 5;}

int g_n_out(void) { return 1;}

const char* g_name_in(int i){
  switch (i) {
  case 0: return "t";
  case 1: return "x";
  case 2: return "u";
  case 3: return "th";
  case 4: return "e";
  default: return 0;
  }
}

const char* g_name_out(int i){
  switch (i) {
  case 0: return "g";
  default: return 0;
  }
}

const int* g_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  case 2: return s0;
  case 3: return s0;
  case 4: return s0;
  default: return 0;
  }
}

const int* g_sparsity_out(int i) {
  switch (i) {
  case 0: return s0;
  default: return 0;
  }
}

int g_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 2;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_g(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  int i, j;
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"g\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"g\" failed. Too many output arguments (%d, max 1)", resc);
  int *iw = 0;
  real_t w[9];
  const real_t* arg[5] = {0};
  if (--argc>=0) arg[0] = from_mex(argv[0], w, s0, w+7);
  if (--argc>=0) arg[1] = from_mex(argv[1], w+1, s1, w+7);
  if (--argc>=0) arg[2] = from_mex(argv[2], w+3, s0, w+7);
  if (--argc>=0) arg[3] = from_mex(argv[3], w+4, s0, w+7);
  if (--argc>=0) arg[4] = from_mex(argv[4], w+5, s0, w+7);
  real_t* res[1] = {0};
  --resc;
  res[0] = w+6;
  i = g(arg, res, iw, w+7, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"g\" failed.");
  if (res[0]) resv[0] = to_mex(s0, res[0]);
}
#endif

/* dgdx */
int dgdx(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=arg[1] ? arg[1][0] : 0;
  a0=(a0+a0);
  real_t a1=2.0000000000000001e-001;
  a1=(a1*a0);
  if (res[0]!=0) res[0][0]=a1;
  return 0;
}

void dgdx_incref(void) {
}

void dgdx_decref(void) {
}

int dgdx_n_in(void) { return 5;}

int dgdx_n_out(void) { return 1;}

const char* dgdx_name_in(int i){
  switch (i) {
  case 0: return "t";
  case 1: return "x";
  case 2: return "u";
  case 3: return "th";
  case 4: return "e";
  default: return 0;
  }
}

const char* dgdx_name_out(int i){
  switch (i) {
  case 0: return "dgdx";
  default: return 0;
  }
}

const int* dgdx_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  case 2: return s0;
  case 3: return s0;
  case 4: return s0;
  default: return 0;
  }
}

const int* dgdx_sparsity_out(int i) {
  switch (i) {
  case 0: return s6;
  default: return 0;
  }
}

int dgdx_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 2;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_dgdx(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  int i, j;
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dgdx\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dgdx\" failed. Too many output arguments (%d, max 1)", resc);
  int *iw = 0;
  real_t w[9];
  const real_t* arg[5] = {0};
  if (--argc>=0) arg[0] = from_mex(argv[0], w, s0, w+7);
  if (--argc>=0) arg[1] = from_mex(argv[1], w+1, s1, w+7);
  if (--argc>=0) arg[2] = from_mex(argv[2], w+3, s0, w+7);
  if (--argc>=0) arg[3] = from_mex(argv[3], w+4, s0, w+7);
  if (--argc>=0) arg[4] = from_mex(argv[4], w+5, s0, w+7);
  real_t* res[1] = {0};
  --resc;
  res[0] = w+6;
  i = dgdx(arg, res, iw, w+7, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dgdx\" failed.");
  if (res[0]) resv[0] = to_mex(s6, res[0]);
}
#endif

/* dgde */
int dgde(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=1.;
  if (res[0]!=0) res[0][0]=a0;
  return 0;
}

void dgde_incref(void) {
}

void dgde_decref(void) {
}

int dgde_n_in(void) { return 5;}

int dgde_n_out(void) { return 1;}

const char* dgde_name_in(int i){
  switch (i) {
  case 0: return "t";
  case 1: return "x";
  case 2: return "u";
  case 3: return "th";
  case 4: return "e";
  default: return 0;
  }
}

const char* dgde_name_out(int i){
  switch (i) {
  case 0: return "dgde";
  default: return 0;
  }
}

const int* dgde_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  case 2: return s0;
  case 3: return s0;
  case 4: return s0;
  default: return 0;
  }
}

const int* dgde_sparsity_out(int i) {
  switch (i) {
  case 0: return s0;
  default: return 0;
  }
}

int dgde_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 1;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_dgde(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  int i, j;
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dgde\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dgde\" failed. Too many output arguments (%d, max 1)", resc);
  int *iw = 0;
  real_t w[9];
  const real_t* arg[5] = {0};
  if (--argc>=0) arg[0] = from_mex(argv[0], w, s0, w+7);
  if (--argc>=0) arg[1] = from_mex(argv[1], w+1, s1, w+7);
  if (--argc>=0) arg[2] = from_mex(argv[2], w+3, s0, w+7);
  if (--argc>=0) arg[3] = from_mex(argv[3], w+4, s0, w+7);
  if (--argc>=0) arg[4] = from_mex(argv[4], w+5, s0, w+7);
  real_t* res[1] = {0};
  --resc;
  res[0] = w+6;
  i = dgde(arg, res, iw, w+7, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dgde\" failed.");
  if (res[0]) resv[0] = to_mex(s0, res[0]);
}
#endif

/* dgdth */
int dgdth(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  return 0;
}

void dgdth_incref(void) {
}

void dgdth_decref(void) {
}

int dgdth_n_in(void) { return 5;}

int dgdth_n_out(void) { return 1;}

const char* dgdth_name_in(int i){
  switch (i) {
  case 0: return "t";
  case 1: return "x";
  case 2: return "u";
  case 3: return "th";
  case 4: return "e";
  default: return 0;
  }
}

const char* dgdth_name_out(int i){
  switch (i) {
  case 0: return "dgdth";
  default: return 0;
  }
}

const int* dgdth_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  case 2: return s0;
  case 3: return s0;
  case 4: return s0;
  default: return 0;
  }
}

const int* dgdth_sparsity_out(int i) {
  switch (i) {
  case 0: return s7;
  default: return 0;
  }
}

int dgdth_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_dgdth(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  int i, j;
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dgdth\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dgdth\" failed. Too many output arguments (%d, max 1)", resc);
  int *iw = 0;
  real_t w[8];
  const real_t* arg[5] = {0};
  if (--argc>=0) arg[0] = from_mex(argv[0], w, s0, w+6);
  if (--argc>=0) arg[1] = from_mex(argv[1], w+1, s1, w+6);
  if (--argc>=0) arg[2] = from_mex(argv[2], w+3, s0, w+6);
  if (--argc>=0) arg[3] = from_mex(argv[3], w+4, s0, w+6);
  if (--argc>=0) arg[4] = from_mex(argv[4], w+5, s0, w+6);
  real_t* res[1] = {0};
  --resc;
  res[0] = w+6;
  i = dgdth(arg, res, iw, w+6, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"dgdth\" failed.");
  if (res[0]) resv[0] = to_mex(s7, res[0]);
}
#endif

/* ddgdxdth */
int ddgdxdth(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  return 0;
}

void ddgdxdth_incref(void) {
}

void ddgdxdth_decref(void) {
}

int ddgdxdth_n_in(void) { return 5;}

int ddgdxdth_n_out(void) { return 1;}

const char* ddgdxdth_name_in(int i){
  switch (i) {
  case 0: return "t";
  case 1: return "x";
  case 2: return "u";
  case 3: return "th";
  case 4: return "e";
  default: return 0;
  }
}

const char* ddgdxdth_name_out(int i){
  switch (i) {
  case 0: return "ddgdxdth";
  default: return 0;
  }
}

const int* ddgdxdth_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  case 2: return s0;
  case 3: return s0;
  case 4: return s0;
  default: return 0;
  }
}

const int* ddgdxdth_sparsity_out(int i) {
  switch (i) {
  case 0: return s8;
  default: return 0;
  }
}

int ddgdxdth_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_ddgdxdth(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  int i, j;
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"ddgdxdth\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"ddgdxdth\" failed. Too many output arguments (%d, max 1)", resc);
  int *iw = 0;
  real_t w[8];
  const real_t* arg[5] = {0};
  if (--argc>=0) arg[0] = from_mex(argv[0], w, s0, w+6);
  if (--argc>=0) arg[1] = from_mex(argv[1], w+1, s1, w+6);
  if (--argc>=0) arg[2] = from_mex(argv[2], w+3, s0, w+6);
  if (--argc>=0) arg[3] = from_mex(argv[3], w+4, s0, w+6);
  if (--argc>=0) arg[4] = from_mex(argv[4], w+5, s0, w+6);
  real_t* res[1] = {0};
  --resc;
  res[0] = w+6;
  i = ddgdxdth(arg, res, iw, w+6, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"ddgdxdth\" failed.");
  if (res[0]) resv[0] = to_mex(s8, res[0]);
}
#endif


#ifdef MATLAB_MEX_FILE
void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  char buf[10];
  int buf_ok = --argc >= 0 && !mxGetString(*argv++, buf, sizeof(buf));
  if (!buf_ok) {
    /* name error */
  } else if (strcmp(buf, "fd")==0) {
    return mex_fd(resc, resv, argc, argv);
  } else if (strcmp(buf, "dfddx")==0) {
    return mex_dfddx(resc, resv, argc, argv);
  } else if (strcmp(buf, "dfddw")==0) {
    return mex_dfddw(resc, resv, argc, argv);
  } else if (strcmp(buf, "dfddth")==0) {
    return mex_dfddth(resc, resv, argc, argv);
  } else if (strcmp(buf, "ddfddxdth")==0) {
    return mex_ddfddxdth(resc, resv, argc, argv);
  } else if (strcmp(buf, "g")==0) {
    return mex_g(resc, resv, argc, argv);
  } else if (strcmp(buf, "dgdx")==0) {
    return mex_dgdx(resc, resv, argc, argv);
  } else if (strcmp(buf, "dgde")==0) {
    return mex_dgde(resc, resv, argc, argv);
  } else if (strcmp(buf, "dgdth")==0) {
    return mex_dgdth(resc, resv, argc, argv);
  } else if (strcmp(buf, "ddgdxdth")==0) {
    return mex_ddgdxdth(resc, resv, argc, argv);
  }
  mexErrMsgTxt("First input should be a command string. Possible values: 'fd' 'dfddx' 'dfddw' 'dfddth' 'ddfddxdth' 'g' 'dgdx' 'dgde' 'dgdth' 'ddgdxdth'");
}
#endif
#ifdef __cplusplus
} /* extern "C" */
#endif
