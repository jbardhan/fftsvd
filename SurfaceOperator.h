/* Typedef */

typedef enum { INSIDE_BIOMOLECULE, INSIDE_STERN, OUTSIDE_STERN } SLIC_regiontype;

typedef struct _SurfaceOperator {
   BEMKernelType kernel;
   real epsilon;
   real kappa;
   unsigned int mynumpanels;
  unsigned int startphi;  // in the global solution vector, x[startphi] = the first entry of the phi for this surface 
  unsigned int startdphi;  // in the global solution vector, x[startdphi] = the first entry of the dphi/dn for this surface
   Tree tree;
   Tree M3;  // optional tree for M3 operator for this region
   unsigned int numchildren;
   struct _SurfaceOperator** children;

   // added by jpb
   Vector resultInternal;
   Vector resultExternal;
   struct _SurfaceOperator* parent;

   Charge charges;
} _SurfaceOperator;

typedef _SurfaceOperator* SurfaceOperator;

/* Constructors and Destructors */

SurfaceOperator SurfaceOperator_allocate();
void SurfaceOperator_free(SurfaceOperator so);
void SurfaceOperator_initStartIndices(SurfaceOperator so, unsigned int *index);
void SurfaceOperator_joinpotentials(SurfaceOperator so, Vector src, Vector phi, Vector dphi, unsigned int scaleDielectric);
void SurfaceOperator_splitglobal(SurfaceOperator so, Vector src);
void SurfaceOperator_joinglobal(SurfaceOperator so, Vector dest);
void SurfaceOperator_topmultiply(SurfaceOperator so, Vector dest, Vector src);
void SurfaceOperator_multiply(SurfaceOperator so, Vector src);
void SurfaceOperator_resultclear(SurfaceOperator so);
void SurfaceOperator_generatePreconditioner(SurfaceOperator so, Preconditioner P);
void SurfaceOperator_generatePreconditioner_diagonal(SurfaceOperator so, Preconditioner P);
void SurfaceOperator_generatePreconditioner_identity(SurfaceOperator so, Preconditioner P);
void SurfaceOperator_makeRHS(SurfaceOperator so, Vector RHS);
void SurfaceOperator_calculateReactionPotentials(SurfaceOperator so, Vector surfacePotentials, 
						 Vector pointPotentials);
void SurfaceOperator_writematlabfile(char* filename, SurfaceOperator so, unsigned int numtotalpanels);
// added by JPB, 12/16/04 to make it easier to integrate with implicit optimization
void SurfaceOperator_makeRHS_fromVector(SurfaceOperator so, Vector RHS,
                                        Vector srcCharges, unsigned int totalnumcharges);
void SurfaceOperator_solve(SurfaceOperator so, Vector sol, Vector RHS);
void SurfaceOperator_collectPotentials_toVector(SurfaceOperator so, Vector phi_r, Vector sol);

