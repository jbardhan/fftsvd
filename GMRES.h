void GMRES_cap(Tree tree, Preconditioner preconditioner, Vector rhs, Vector sol, real tol);
void GMRES_solv(Tree tree, Preconditioner preconditioner, Vector rhs, Vector sol, real tol, real idiel, real odiel);
void GMRES_solv_ecf(Tree tree, Preconditioner preconditioner, Vector rhs, Vector sol, real tol, real idiel, real odiel);
void GMRES_solv_ecf_qual(Tree tree, Preconditioner preconditioner, Vector rhs, Vector sol, real tol, real idiel, real odiel);
void GMRES_SurfaceOperator(SurfaceOperator so, Preconditioner preconditioner, Vector rhs, Vector sol, unsigned int numpanels, real tol);
void GMRES_QualocationOperator(QualocationOperator qo, Preconditioner preconditioner, Vector rhs, Vector sol, unsigned int numpanels, real tol);

void Solv_ecf_multiply_qual(Vector b, Tree tree, Vector x, real idiel, real odiel);
void Solv_multiply(Vector b, Tree tree, Vector x, real idiel, real odiel);

extern unsigned int num_GMRES_iter;

