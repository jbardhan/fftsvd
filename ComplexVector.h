/* Typedef */

typedef complexreal* ComplexVector;

/* Constructors and Destructors */

ComplexVector ComplexVector_allocate(unsigned int length);
void ComplexVector_free(ComplexVector vector);

/* Initialization and Copying */

void ComplexVector_copy(ComplexVector vectordest, ComplexVector vectorsrc, unsigned int length);
void ComplexVector_zero(ComplexVector vector, unsigned int length);

/* Arithmetic Operations */

void ComplexVector_addelementmultiplyvector(ComplexVector vector1, ComplexVector vector2, ComplexVector vector3, unsigned int length);

void ComplexVector_writefile(char* filename, ComplexVector vector, unsigned int length);
