#ifndef DSEUPD_HPP
#define DSEUPD_HPP

extern "C"
{
  void dseupd_(bool *rvec, char *HowMny, bool *select,
               double *d, double *Z, int *ldz,
               double *sigma, char *bmat, int *n,
               char *which, int *nev, double *tol,
               double *resid, int *ncv, double *V,
               int *ldv, int *iparam, int *ipntr,
               double *workd, double *workl,
               int *lworkl, int *info);
}

void dseupd(bool &rvec, char &HowMny, bool *select,
            double *d, double *Z, int &ldz,
            double &sigma, char &bmat, int &n,
            char *which, int &nev, double &tol,
            double *resid, int &ncv, double *V,
            int &ldv, int *iparam, int *ipntr,
            double *workd, double *workl,
            int &lworkl, int &info);

inline
void dseupd(bool &rvec, char &howmny, bool *select,
            double *d, double *z, int &ldz,
            double &sigma, char &bmat, int &n,
            char *which, int &nev, double &tol,
            double *resid, int &ncv, double *v,
            int &ldv, int *iparam, int *ipntr,
            double *workd, double *workl,
            int &lworkl, int &info)
{
  dseupd_( &rvec, &howmny, select, d, z, &ldz, &sigma, 
           &bmat, &n, which, &nev, &tol, resid, &ncv, 
           v, &ldv, iparam, ipntr, workd, workl, &lworkl, 
           &info );
}

#endif
