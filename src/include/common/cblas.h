// Head file for BLAS


#ifndef __CBLAS_H__
#define __CBLAS_H__


#undef BLAS_INTEGER
#define BLAS_INTEGER int
#undef BLAS_REAL
#define BLAS_REAL float
#undef BLAS_DOUBLEREAL
#define BLAS_DOUBLEREAL double
#undef BLAS_COMPLEX
#define BLAS_COMPLEX void
#undef BLAS_DOUBLECOMPLEX
#define BLAS_DOUBLECOMPLEX void
#undef BLAS_LOGICAL
#define BLAS_LOGICAL int
#undef BLAS_L_FP
#define BLAS_L_FP int*
#undef BLAS_FTNLEN
#define BLAS_FTNLEN int*


/* Subroutine */ void dgemm_(char *transa, char *transb, BLAS_INTEGER *m, BLAS_INTEGER *n, BLAS_INTEGER *k,  BLAS_DOUBLEREAL *alpha,
			     BLAS_DOUBLEREAL *a, BLAS_INTEGER *tda,  BLAS_DOUBLEREAL *b, BLAS_INTEGER *tdb,
			     BLAS_DOUBLEREAL *beta,  BLAS_DOUBLEREAL *c, BLAS_INTEGER *tdc);


#endif /* !__CBLAS_H__ */

