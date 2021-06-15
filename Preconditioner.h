/* Typedef */

typedef struct PreconditionerElement {
   unsigned int row;
   real value;
} PreconditionerElement;

typedef struct _Preconditioner {
   unsigned int rows, columns, nnz;
   PreconditionerElement** P;
   unsigned int* numelements;
   SuperMatrix L;
   SuperMatrix U;
   int* perm_r;
   int* perm_c;   
   SuperLUStat_t stat;
} _Preconditioner;

typedef _Preconditioner* Preconditioner;

/* Constructors and Destructors */

Preconditioner Preconditioner_allocate(unsigned int rows, unsigned int columns);
void Preconditioner_free(Preconditioner preconditioner);

/* Operations */

void Preconditioner_set(Preconditioner preconditioner, unsigned int row, unsigned int column, real value);
void Preconditioner_fill_identity(Preconditioner preconditioner);
void Preconditioner_fill_diagonal_cap(Preconditioner preconditioner, Tree tree);
void Preconditioner_fill_diagonal_solv_ecf_qual(Preconditioner preconditioner, Tree tree, real idiel, real odiel);
void Preconditioner_fill_diagonal_solv_born(Preconditioner preconditioner, Tree tree, real idiel, real odiel);
void Preconditioner_fill_diagonal_solv_RUB(Preconditioner preconditioner, Tree tree, real idiel, real odiel);
void Preconditioner_fill_diagonal_solv_juffersimple(Preconditioner preconditioner, Tree tree, real idiel, real odiel);
void Preconditioner_fill_blockdiagonal_solv(Preconditioner preconditioner, Tree tree, real idiel, real odiel);
void Preconditioner_fill_blockdiagonal_solv_cavity(Preconditioner preconditioner, Panel* panels1, unsigned int numpanels1, Panel* panels2, unsigned int numpanels2, BEMKernelType kerneltype, void* parameters);
void Preconditioner_factor(Preconditioner preconditioner);
void Preconditioner_solve(Vector x, Preconditioner preconditioner, Vector b);
unsigned int Preconditioner_memory(Preconditioner preconditioner);
void Preconditioner_writematlabfile(char* filename, Preconditioner preconditioner);
