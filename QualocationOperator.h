/* Typedef */

typedef struct _QualocationOperator {
   Tree tree;
   Tree M1M3;  // optional combined tree for M1/M3

   real innerdielectric, outerdielectric;

   Charge charges;
} _QualocationOperator;

typedef _QualocationOperator* QualocationOperator;

/* Constructors and Destructors */

QualocationOperator QualocationOperator_allocate();
void QualocationOperator_free(QualocationOperator qo);

void QualocationOperator_multiply(QualocationOperator so, Vector dest, Vector src);
void QualocationOperator_generatePreconditioner(QualocationOperator qo, Preconditioner P);
void QualocationOperator_makeRHS(QualocationOperator qo, Vector RHS);
void QualocationOperator_calculateReactionPotentials(QualocationOperator qo, Vector surfacePotentials, 
						 Vector pointPotentials);
void QualocationOperator_writematlabfile(char* filename, QualocationOperator qo, unsigned int numtotalpanels);
void QualocationOperator_makeRHS_fromVector(QualocationOperator qo, Vector RHS,
                                        Vector srcCharges, unsigned int totalnumcharges);
void QualocationOperator_solve(QualocationOperator qo, Vector sol, Vector RHS);
void QualocationOperator_collectPotentials_toVector(QualocationOperator qo, Vector phi_r, Vector sol);


