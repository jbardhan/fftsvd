#include "FFTSVD.h"
 
// macro to get array reference from cube indices
// the indices are in the padded cube (frame of reference)
#define CREF(N, i, j, k) ((k) + N*((j) + N* (i)))

// n is the number of grid points we're using, ie unpadded
// shift is the distance vector
// gridspacing, well, duh
void FFT_calcDiagonalTranslationOperator(unsigned int level, unsigned int n, Vector3D shift, real gridspacing, BEMKernelType kerneltype, void* parameters, ComplexSVector translate) {
   unsigned int N = 2*n-1;
   SVector data = SVector_allocate(N * N * N);
   int i, j, k, li, lj, lk;
   Vector3D zero = Vector3D_allocate();
   Vector3D point = Vector3D_allocate();

   for (i=0; i < N; i++) {
      li = -i; if (i > N/2) li = N - i;
      for (j=0; j < N; j++) {
         lj = -j; if (j > N/2) lj = N - j;
         for (k=0; k < N; k++) {
            lk = -k; if (k > N/2) lk = N - k;
            point->x = gridspacing * (shift->x - li);
            point->y = gridspacing * (shift->y - lj);
            point->z = gridspacing * (shift->z - lk);
            data[CREF(N, i, j, k)] = GreensFunction(point, zero, kerneltype, parameters);
         }
      }
   }

   Vector3D_free(zero);
   Vector3D_free(point);

#ifdef SREAL_IS_DOUBLE
   {
      static fftw_plan plans[10] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
#ifdef OMP
#pragma omp critical (fftw)
#endif
      if (plans[n] == NULL) {
         plans[n] = fftw_plan_dft_r2c_3d(N, N, N, data, translate, FFTW_ESTIMATE | FFTW_UNALIGNED);
      }
      fftw_execute_dft_r2c(plans[n], data, translate);
      //fftw_destroy_plan(plan);
   }
#else
#ifdef SREAL_IS_FLOAT
   {
      static fftwf_plan plans[10] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
#ifdef OMP
#pragma omp critical (fftw)
#endif
      if (plans[n] == NULL)
         plans[n] = fftwf_plan_dft_r2c_3d(N, N, N, data, translate, FFTW_ESTIMATE | FFTW_UNALIGNED);
      fftwf_execute_dft_r2c(plans[n], data, translate);
      //fftwf_destroy_plan(plan);
   }
#endif
#endif

   SVector_free(data);
}
 
// this function converts our UNPADDED grid charges (array size n^3)
// into PADDED grid Fourier coefficients(array size (2*n-1)^3)
void FFT_forwardGridTransform(unsigned int level, unsigned int n, SVector gridCharges, ComplexSVector paddedFourierCharges) {
   unsigned int N = 2*n-1;
   unsigned int offset = n/2;
   SVector paddedGridCharges = SVector_allocate(N * N * N);
   int i, j, k;

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         for (k = 0; k < n; k++) {
            paddedGridCharges[CREF(N,(i+offset),(j+offset),(k+offset))] =
               gridCharges[CREF(n, i, j, k)];
         }
      }
   }

#ifdef SREAL_IS_DOUBLE
   {
      static fftw_plan plans[10] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
#ifdef OMP
#pragma omp critical (fftw)
#endif
      if (plans[n] == NULL)
         plans[n] = fftw_plan_dft_r2c_3d(N, N, N, paddedGridCharges, paddedFourierCharges, FFTW_ESTIMATE | FFTW_UNALIGNED);
      fftw_execute_dft_r2c(plans[n], paddedGridCharges, paddedFourierCharges);
      //fftw_destroy_plan(plan);
   }
#else
#ifdef SREAL_IS_FLOAT
   {
      static fftwf_plan plans[10] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
#ifdef OMP
#pragma omp critical (fftw)
#endif
      if (plans[n] == NULL)
         plans[n] = fftwf_plan_dft_r2c_3d(N, N, N, paddedGridCharges, paddedFourierCharges, FFTW_ESTIMATE | FFTW_UNALIGNED);
      fftwf_execute_dft_r2c(plans[n], paddedGridCharges, paddedFourierCharges);
      //fftwf_destroy_plan(plan);
   }
#endif
#endif 

   SVector_free(paddedGridCharges);
}
 
// this function takes our PADDED grid Fourier coefficients (array of size (2*n-1)^3)
// and converts it back into grid potential space
void FFT_backwardGridTransform(unsigned int level, unsigned int n, ComplexSVector paddedFourierResponse, SVector gridPotentials) {
   unsigned int N = 2*n-1;
   unsigned int offset = n/2;
   SVector paddedGridPotentials = SVector_allocate(N * N * N);
   size_t i, j, k;

#ifdef SREAL_IS_DOUBLE
   static fftw_plan plans[10] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
#ifdef OMP
#pragma omp critical (fftw)
#endif
   if (plans[n] == NULL)
      plans[n] = fftw_plan_dft_c2r_3d(N, N, N, paddedFourierResponse, paddedGridPotentials, FFTW_ESTIMATE | FFTW_UNALIGNED);
   fftw_execute_dft_c2r(plans[n], paddedFourierResponse, paddedGridPotentials);
   //fftw_destroy_plan(plan);
#else
#ifdef SREAL_IS_FLOAT
#endif
   static fftwf_plan plans[10] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
#ifdef OMP
#pragma omp critical (fftw)
#endif
   if (plans[n] == NULL)
      plans[n] = fftwf_plan_dft_c2r_3d(N, N, N, paddedFourierResponse, paddedGridPotentials, FFTW_ESTIMATE | FFTW_UNALIGNED);
   fftwf_execute_dft_c2r(plans[n], paddedFourierResponse, paddedGridPotentials);
   //fftwf_destroy_plan(plan);
#endif 

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         for (k = 0; k < n; k++) {
            gridPotentials[CREF(n, i, j, k)] =
               paddedGridPotentials[CREF(N, (i+offset), (j+offset), (k+offset))] / (N * N * N);
         }
      }
   }
 
   SVector_free(paddedGridPotentials);
}
