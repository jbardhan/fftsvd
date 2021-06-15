/* Typedef */

typedef real** Matrix;

/* Constructors and Destructors */

Matrix Matrix_allocate(unsigned int rows, unsigned int columns);
void Matrix_free(Matrix matrix);

/* Initialization and Copying */

void Matrix_copy(Matrix matrixdest, Matrix matrixsrc, unsigned int rows, unsigned int columns);
void Matrix_copypiece(Matrix matrixdest, unsigned int destStartRow, unsigned int destStartCol,
                      Matrix matrixsrc, unsigned int srcStartRow, unsigned int srcStartCol,
                      unsigned int numrows, unsigned int numcols);
void Matrix_zero(Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_identity(Matrix matrix, unsigned int rows, unsigned int columns);

void Matrix_setColumnFromVector(Matrix matrix, unsigned int index, Vector column, unsigned int rows, unsigned int columns);
void Matrix_getColumnIntoVector(Matrix matrix, unsigned int index, Vector column, unsigned int rows, unsigned int columns);

/* Linear Algebra Operations */

void Matrix_scale(Matrix A, real factor, unsigned int rows, unsigned int columns);
void Matrix_add(Matrix B, Matrix A, Matrix X, unsigned int rows, unsigned int columns);
void Matrix_stats(Matrix matrix, unsigned int rows, unsigned int columns);
real Matrix_norm(Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_multiplyvector(Vector b, Matrix A, Vector x, unsigned int rows, unsigned int columns);
void Matrix_multiplyvector_transpose(Vector b, Matrix A, Vector x, unsigned int rows, unsigned int columns);
void Matrix_multiplymatrix(Matrix B, Matrix A, Matrix X, unsigned int rowsA, unsigned int columnsA, unsigned int columnsX);
void Matrix_multiplymatrix_transpose(Matrix B, Matrix A, Matrix X, unsigned int rowsA, unsigned int columnsA, unsigned int columnsX);
void Matrix_multiplytranspose_matrix(Matrix B, Matrix A, Matrix X, unsigned int rowsA, unsigned int columnsA, unsigned int rowsX);
void Matrix_rowbasis_pivotedmgs(Matrix* U, unsigned int* rowrank, Matrix matrix, unsigned int rows, unsigned int columns, real epsilon);
void Matrix_columnbasis_pivotedmgs(Matrix* V, unsigned int* columnrank, Matrix matrix, unsigned int rows, unsigned int columns, real epsilon);
void Matrix_transpose(Matrix* matrix, unsigned int rows, unsigned int columns);
void Matrix_diff(Matrix X, Matrix A, Matrix B, unsigned int rows, unsigned int columns);
void Matrix_pseudoinverse_droptol(Matrix XI, Matrix X, unsigned int rows, unsigned int columns, real droptol);
void Matrix_pseudoinverse(Matrix XI, Matrix X, unsigned int rows, unsigned int columns);
void Matrix_lsq_solve(Vector x, Matrix A, Vector b, unsigned int rows, unsigned int columns);
void Matrix_columnbasis_check(Matrix V, unsigned int rows, unsigned int columns);
void Matrix_eigendecomposition(Matrix A, Matrix V, Vector d, unsigned int rows, unsigned int columns);
//       Matrix_eigendecomposition symmetrizes A first (in a temporary variable)!!
void Matrix_QR(Matrix A, Matrix Q, Matrix R, unsigned int rows, unsigned int columns);

/* I/O */

void Matrix_writefile(char* filename, Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_readfile(Matrix matrix, FILE* file, unsigned int rows, unsigned int columns);
void Matrix_print(Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_writebinary(FILE* file, Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_readbinary(Matrix matrix, FILE* file, unsigned int rows, unsigned int columns);
