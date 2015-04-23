#ifndef SIMOL_DSAUPD_HPP
#define SIMOL_DSAUPD_HPP

namespace simol
{

  extern "C"
  {
    void dsaupd_(int *ido, char *bmat, int *n, char *which,
                 int *nev, double *tol, double *resid,
                 int *ncv, double *v, int *ldv,
                 int *iparam, int *ipntr, double *workd,
                 double *workl, int *lworkl, int *info);
  }

  void dsaupd(int &ido, char &bmat, int &n, char *which,
              int &nev, double &tol, double *resid,
              int &ncv, double *v, int &ldv,
              int *iparam, int *ipntr, double *workd,
              double *workl, int &lworkl, int &info);

  inline
  void dsaupd(int &ido, char &bmat, int &n, char *which,
              int &nev, double &tol, double *resid,
              int &ncv, double *v, int &ldv,
              int *iparam, int *ipntr, double *workd,
              double *workl, int &lworkl, int &info)
  {
    dsaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, 
            &ncv, v, &ldv, iparam, ipntr, workd, workl, 
           &lworkl, &info);
  }
}

#endif
