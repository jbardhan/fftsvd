#include "FFTSVD.h"

/* Constructors and Destructors */

SVector SVector_allocate(unsigned int length) {
   return (SVector)calloc(length, sizeof(sreal));
}

void SVector_free(SVector vector) {
   free(vector);
}

/* Initialization and Copying */

void SVector_copy(SVector vectordest, SVector vectorsrc, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vectordest[i] = vectorsrc[i];
}

void SVector_copypiece(SVector vectordest, unsigned int destStart, SVector vectorsrc, unsigned int srcStart,unsigned int length) {
   unsigned int i;
#pragma ivdep
   for (i = 0; i < length; i++)
      vectordest[destStart + i] = vectorsrc[srcStart + i];
}

void SVector_copyscaledpiece(SVector vectordest, unsigned int destStart, SVector vectorsrc, unsigned int srcStart,unsigned int length, sreal scalefactor) {
   unsigned int i;
#pragma ivdep
   for (i = 0; i < length; i++)
      vectordest[destStart + i] = scalefactor * vectorsrc[srcStart + i];
}

void SVector_zero(SVector vector, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector[i] = (sreal)0;
}

/* Arithmetic Operations */

void SVector_scale(SVector vector, sreal scale, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector[i] *= scale;
}

void SVector_subtractscaledvector(SVector vector1, sreal scale, SVector vector2, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector1[i] -= scale * vector2[i];
}

void SVector_addscaledvector(SVector vector1, sreal scale, SVector vector2, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector1[i] += scale * vector2[i];
}

void SVector_addvector(SVector vector1, SVector vector2, unsigned int length) {
   unsigned int i;
 
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector1[i] += vector2[i];
}

void SVector_subtractvector(SVector vector1, SVector vector2, unsigned int length) {
   unsigned int i;
 
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector1[i] -= vector2[i];
}

/* Linear Algebra Operations */

sreal SVector_dot(SVector vector1, SVector vector2, unsigned int length) {
   unsigned int i;
   sreal dot = (sreal)0;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      dot += vector1[i] * vector2[i];
      
   return dot;
}

sreal SVector_norm(SVector vector, unsigned int length) {
   unsigned int i;
   sreal norm = (sreal)0;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      norm += vector[i] * vector[i];
   
   norm = sqrt(norm);
   
   return norm;
}

void SVector_normalize(SVector vector, unsigned int length) {
   SVector_scale(vector, (sreal)1 / SVector_norm(vector, length), length);
}
