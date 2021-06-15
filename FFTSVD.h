#ifndef FFTSVD_H
#define FFTSVD_H

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <sys/stat.h>
#ifndef _AIX
#include <complex.h>
#endif

#ifndef _MALLOC
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

#include "fftw3.h"

#define PRINT_GMRES_RESIDUALS

#define REAL_IS_DOUBLE  /* Size of real for integration */
#undef REAL_IS_FLOAT
#undef SREAL_IS_DOUBLE  /* Size of real for storage */
#define SREAL_IS_FLOAT

/* SuperLU headers */

#ifdef REAL_IS_FLOAT
#include "slu_sdefs.h"
#else
#ifdef REAL_IS_DOUBLE
#include "slu_ddefs.h"
#endif
#endif

/* Our AIX machine is not quite C99 compliant, this helps */

#ifdef _AIX
#define complex _Complex
#define creal __real__
#define cimag __imag__
#endif

#ifdef REAL_IS_DOUBLE
typedef double real;
typedef double complex complexreal;
#else
#ifdef REAL_IS_FLOAT
typedef float real;
typedef float complex complexreal;
#endif
#endif

#ifdef SREAL_IS_DOUBLE
typedef double sreal;
typedef double complex complexsreal;
#else
#ifdef SREAL_IS_FLOAT
typedef float sreal;
typedef float complex complexsreal;
#endif
#endif

#define LOCAL 1  /* Only adjacent cubes are nearest neighbors, as in 
                    Greengard or Nabors, set to 2 to have next-nearest 
                    neighbors on the direct list */
#define MAX_QUADRATURE_ORDER 1024  /* Remove */
#undef ACCELERATED_SAMPLING  /* Accelerate sampling for dominant 
                                sources/responses using the FFT, not 
                                recommended for high accuracy */
#undef ACCELERATED_K  /* Compute dominant source to dominant response 
                         "translation" matrices using the FFT rather 
                         than using direct integrations.  Only 
                          meaningful if ADAPTIVE is on */
#undef ADAPTIVE  /* Use "K" matrices rather than the FFT to represent 
                    dominant source to dominant response "translations".  
                    Improves MV product time dramatically at the cost 
                    of high memory usage. */
#undef POLYNOMIAL  /* Use polynomials for computing projection 
                      and interpolation matrices rather than equivalent 
                      density, just like in pFFT++.  Allows for true
                      kernel-independence. */
#undef SERIALIZE  /* Store compressed matrices on disk instead of in 
                     RAM.  The code will use no memory, but you better 
                     have a fast disk! */

typedef enum { CONSTANT_KERNEL, X_KERNEL, Y_KERNEL, Z_KERNEL, POISSON_KERNEL, HELMHOLTZ_KERNEL, DESINGULARIZED_HELMHOLTZ_KERNEL, LJ_KERNEL, LJ12_KERNEL, LJ6_KERNEL, MONOMIAL_KERNEL, INVERSEPOWER_KERNEL, GHOSH_KERNEL, GRYCUK_KERNEL, LESYNG_KERNEL, VOLUME_KERNEL, XGB_KERNEL, JUFFER_KERNEL, BIBEE_P_KERNEL, BIBEE_CFA_KERNEL} BEMKernelType;
typedef enum { SINGLE_LAYER_INT, DOUBLE_LAYER_INT, SINGLE_AND_DOUBLE_LAYER_INT, NORMDERIV_SINGLE_LAYER_INT, NORMDERIV_DOUBLE_LAYER_INT } BEMLayerType;

#undef ENABLE_GALERKIN
#define GALERKIN_ORDER 1
typedef enum { COLLOCATION, QUALOCATION, BASIS_INNER, TEST_INNER, POINT_CHARGE, E_FIELD_TEST_ONLY} BEMDiscretizationType;
// rules 1, 3, 5, 10, 13, 18 

#undef SCATTER
#ifdef SCATTER
#define _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
extern unsigned int resolution;
#else
#undef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
#endif

#include "Utility.h"
#include "Vector.h"
#include "ComplexVector.h"
#include "Matrix.h"
#include "SVector.h"
#include "ComplexSVector.h"
#include "SMatrix.h"
#include "Vector3D.h"

#include "FlatPanel.h"
#include "GST.h"
#include "TOR.h"
#include "Panel.h"

#include "QuadratureRule.h"
#include "Integration.h"
#include "VertFace.h"
#include "QUI.h"
#include "Charge.h"
#include "LJparameters.h"
#include "Cube.h"
#include "EquivDensity.h"
#include "GreensFunction.h"
#include "Tree.h"
#include "Preconditioner.h"
#include "SurfaceOperator.h"
#include "QualocationOperator.h"

#define MAXITERATIONS 1000

#include "GMRES.h"
#include "FFT.h"
#ifdef POLYNOMIAL
#include "Polynomial.h"
#endif

#include "FlatIntegration.h"
#include "calcpc_GST.h"
#include "calcpc_TOR.h"

#endif
