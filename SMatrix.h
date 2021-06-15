/* Typedef */

typedef sreal** SMatrix;

/* Constructors and Destructors */

SMatrix SMatrix_allocate(unsigned int rows, unsigned int columns);
void SMatrix_free(SMatrix matrix);

/* Initialization and Copying */

void SMatrix_copy(SMatrix matrixdest, SMatrix matrixsrc, unsigned int rows, unsigned int columns);
void SMatrix_copypiece(SMatrix matrixdest, unsigned int destStartRow, unsigned int destStartCol,
                      SMatrix matrixsrc, unsigned int srcStartRow, unsigned int srcStartCol,
                      unsigned int numrows, unsigned int numcols);
void SMatrix_zero(SMatrix matrix, unsigned int rows, unsigned int columns);
void SMatrix_identity(SMatrix matrix, unsigned int rows, unsigned int columns);
void SMatrix_Matrix(SMatrix smatrix, Matrix matrix, unsigned int rows, unsigned int columns);

/* Linear Algebra Operations */

void SMatrix_stats(SMatrix matrix, unsigned int rows, unsigned int columns);
sreal SMatrix_norm(SMatrix matrix, unsigned int rows, unsigned int columns);
void SMatrix_multiplyvector(SVector b, SMatrix A, SVector x, unsigned int rows, unsigned int columns);
void SMatrix_multiplyvector_transpose(SVector b, SMatrix A, SVector x, unsigned int rows, unsigned int columns);
void SMatrix_multiplymatrix(SMatrix B, SMatrix A, SMatrix X, unsigned int rowsA, unsigned int columnsA, unsigned int columnsX);
void SMatrix_multiplymatrix_transpose(SMatrix B, SMatrix A, SMatrix X, unsigned int rowsA, unsigned int columnsA, unsigned int columnsX);
void SMatrix_multiplytranspose_matrix(SMatrix B, SMatrix A, SMatrix X, unsigned int rowsA, unsigned int columnsA, unsigned int rowsX);
void SMatrix_rowbasis_pivotedmgs(SMatrix* U, unsigned int* rowrank, SMatrix matrix, unsigned int rows, unsigned int columns, sreal epsilon);
void SMatrix_columnbasis_pivotedmgs(SMatrix* V, unsigned int* columnrank, SMatrix matrix, unsigned int rows, unsigned int columns, sreal epsilon);
void SMatrix_transpose(SMatrix* matrix, unsigned int rows, unsigned int columns);
void SMatrix_diff(SMatrix X, SMatrix A, SMatrix B, unsigned int rows, unsigned int columns);
void SMatrix_pseudoinverse_droptol(SMatrix XI, SMatrix X, unsigned int rows, unsigned int columns, sreal droptol);
void SMatrix_pseudoinverse(SMatrix XI, SMatrix X, unsigned int rows, unsigned int columns);
void SMatrix_lsq_solve(SVector x, SMatrix A, SVector b, unsigned int rows, unsigned int columns);
void SMatrix_columnbasis_check(SMatrix V, unsigned int rows, unsigned int columns);

/* I/O */

void SMatrix_writefile(char* filename, SMatrix matrix, unsigned int rows, unsigned int columns);
void SMatrix_print(SMatrix matrix, unsigned int rows, unsigned int columns);
void SMatrix_writebinary(FILE* file, SMatrix matrix, unsigned int rows, unsigned int columns);
void SMatrix_readbinary(SMatrix matrix, FILE* file, unsigned int rows, unsigned int columns);
