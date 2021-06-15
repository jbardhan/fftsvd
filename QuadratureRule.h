/* Typedef */

typedef struct {
  /* basic integration */
  unsigned int order;
  real *x;
  real *w;
                                                                                
  /* 7 point quadrature */
  real dNdKsi[6][7];
  real dNdEta[6][7];
  real SPW[7];
} _QuadratureRule;

typedef _QuadratureRule* QuadratureRule;

/* Constructors and Destructors */

QuadratureRule QuadratureRule_allocate(unsigned int order);
void QuadratureRule_free(QuadratureRule qr);

/* Operations */

void QuadratureRule_generatepoints(QuadratureRule qr, real lower, real upper);
