/* Typedef */

typedef real* Vector;

/* Constructors and Destructors */

Vector Vector_allocate(unsigned int length);
void Vector_free(Vector vector);

/* Initialization and Copying */

void Vector_copy(Vector vectordest, Vector vectorsrc, unsigned int length);
void Vector_copypiece(Vector vectordest, unsigned int destStart, Vector vectorsrc, unsigned int srcStart, unsigned int length);
void Vector_copyscaledpiece(Vector vectordest, unsigned int destStart, Vector vectorsrc, unsigned int srcStart, unsigned int length, real scalefactor);
void Vector_zero(Vector vector, unsigned int length);

/* Arithmetic Operations */

void Vector_scale(Vector vector, real scale, unsigned int length);
void Vector_subtractscaledvector(Vector vector1, real scale, Vector vector2, unsigned int length);
void Vector_addscaledvector(Vector vector1, real scale, Vector vector2, unsigned int length);
void Vector_addvector(Vector vector1, Vector vector2, unsigned int length);
void Vector_subtractvector(Vector vector1, Vector vector2, unsigned int length);
void Vector_pointwisescale(Vector dest, Vector scalefactors, unsigned int length);

/* Linear Algebra Operations */

real Vector_dot(Vector vector1, Vector vector2, unsigned int length);
real Vector_norm(Vector vector, unsigned int length);
void Vector_normalize(Vector vector, unsigned int length);

/* Input/Output Operations */
void Vector_writefile(char* filename, Vector vector, unsigned int length);

