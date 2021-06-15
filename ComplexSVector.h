/* Typedef */

typedef complexsreal* ComplexSVector;

/* Constructors and Destructors */

ComplexSVector ComplexSVector_allocate(unsigned int length);
void ComplexSVector_free(ComplexSVector vector);

/* Initialization and Copying */

void ComplexSVector_copy(ComplexSVector vectordest, ComplexSVector vectorsrc, unsigned int length);
void ComplexSVector_zero(ComplexSVector vector, unsigned int length);

/* Arithmetic Operations */

void ComplexSVector_addelementmultiplyvector(ComplexSVector vector1, ComplexSVector vector2, ComplexSVector vector3, unsigned int length);

void ComplexSVector_writefile(char* filename, ComplexSVector vector, unsigned int length);
