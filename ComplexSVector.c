#include "FFTSVD.h"
#ifdef INTEL_C
#include "ipps.h"
#endif

/* Constructors and Destructors */

ComplexSVector ComplexSVector_allocate(unsigned int length) {
#ifdef INTEL_C
#ifdef SREAL_IS_DOUBLE
   ComplexSVector vc = (ComplexSVector)ippsMalloc_64fc(length);
#else
   ComplexSVector vc = (ComplexSVector)ippsMalloc_32fc(length);
#endif
   memset(vc, 0, length * sizeof(complexsreal));
   return vc;
#else
   return (ComplexSVector)calloc(length, sizeof(complexsreal));
#endif
}

void ComplexSVector_free(ComplexSVector vector) {
#ifdef INTEL_C
   ippsFree(vector);
#else
   free(vector);
#endif
}

/* Initialization and Copying */

void ComplexSVector_copy(ComplexSVector vectordest, ComplexSVector vectorsrc, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vectordest[i] = vectorsrc[i];
}

void ComplexSVector_zero(ComplexSVector vector, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++) {
      vector[i] = (complexsreal)0;
   }
}

/* Arithmetic Operations */

void ComplexSVector_addelementmultiplyvector(ComplexSVector vector1, ComplexSVector vector2, ComplexSVector vector3, unsigned int length) {
   unsigned int i;

#ifdef INTEL_C
#ifdef SREAL_IS_FLOAT
   ippsAddProduct_32fc((Ipp32fc*)vector2, (Ipp32fc*)vector3, (Ipp32fc*)vector1, length);
#else
#ifdef SREAL_IS_DOUBLE
   ippsAddProduct_64fc((Ipp64fc*)vector2, (Ipp64fc*)vector3, (Ipp64fc*)vector1, length);
#endif
#endif
#else
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector1[i] += vector2[i] * vector3[i];
#endif
}

void ComplexSVector_writefile(char* filename, ComplexSVector vector, unsigned int length) {
   unsigned int i;
   FILE* file = NULL;

   file = fopen(filename, "w");

#ifndef _AIX
   for (i = 0; i < length; i++)
      fprintf(file, "%20.15e %20.15e\n", creal(vector[i]), cimag(vector[i]));
#endif

   fclose(file);
}
