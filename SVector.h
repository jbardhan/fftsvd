/* Typedef */

typedef sreal* SVector;

/* Constructors and Destructors */

SVector SVector_allocate(unsigned int length);
void SVector_free(SVector vector);

/* Initialization and Copying */

void SVector_copy(SVector vectordest, SVector vectorsrc, unsigned int length);
void SVector_copypiece(SVector vectordest, unsigned int destStart, SVector vectorsrc, unsigned int srcStart, unsigned int length);
void SVector_copyscaledpiece(SVector vectordest, unsigned int destStart, SVector vectorsrc, unsigned int srcStart, unsigned int length, sreal scalefactor);
void SVector_zero(SVector vector, unsigned int length);

/* Arithmetic Operations */

void SVector_scale(SVector vector, sreal scale, unsigned int length);
void SVector_subtractscaledvector(SVector vector1, sreal scale, SVector vector2, unsigned int length);
void SVector_addscaledvector(SVector vector1, sreal scale, SVector vector2, unsigned int length);
void SVector_addvector(SVector vector1, SVector vector2, unsigned int length);
void SVector_subtractvector(SVector vector1, SVector vector2, unsigned int length);

/* Linear Algebra Operations */

sreal SVector_dot(SVector vector1, SVector vector2, unsigned int length);
sreal SVector_norm(SVector vector, unsigned int length);
void SVector_normalize(SVector vector, unsigned int length);

