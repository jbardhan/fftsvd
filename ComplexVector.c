#include "FFTSVD.h"
#ifdef INTEL_C
#include "ipps.h"
#endif

/* Constructors and Destructors */

ComplexVector ComplexVector_allocate(unsigned int length) {
#ifdef INTEL_C
#ifdef REAL_IS_DOUBLE
   ComplexVector vc = (ComplexVector)ippsMalloc_64fc(length);
#else
   ComplexVector vc = (ComplexVector)ippsMalloc_32fc(length);
#endif
   memset(vc, 0, length * sizeof(complexreal));
   return vc;
#else
   return (ComplexVector)calloc(length, sizeof(complexreal));
#endif
}

void ComplexVector_free(ComplexVector vector) {
#ifdef INTEL_C
   ippsFree(vector);
#else
   free(vector);
#endif
}

/* Initialization and Copying */

void ComplexVector_copy(ComplexVector vectordest, ComplexVector vectorsrc, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vectordest[i] = vectorsrc[i];
}

void ComplexVector_zero(ComplexVector vector, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++) {
      vector[i] = (complexreal)0;
   }
}

/* Arithmetic Operations */

void ComplexVector_addelementmultiplyvector(ComplexVector vector1, ComplexVector vector2, ComplexVector vector3, unsigned int length) {
   unsigned int i;

#ifdef INTEL_C
#ifdef REAL_IS_FLOAT
   ippsAddProduct_32fc((Ipp32fc*)vector2, (Ipp32fc*)vector3, (Ipp32fc*)vector1, length);
#else
#ifdef REAL_IS_DOUBLE
   ippsAddProduct_64fc((Ipp64fc*)vector2, (Ipp64fc*)vector3, (Ipp64fc*)vector1, length);
#endif
#endif
#else
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector1[i] += vector2[i] * vector3[i];
#endif
}

void ComplexVector_writefile(char* filename, ComplexVector vector, unsigned int length) {
   unsigned int i;
   FILE* file = NULL;

   file = fopen(filename, "w");

#ifndef _AIX
   for (i = 0; i < length; i++)
      fprintf(file, "%20.15e %20.15e\n", creal(vector[i]), cimag(vector[i]));
#endif

   fclose(file);
}
