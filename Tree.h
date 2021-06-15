/* Typedef */

typedef struct _Tree {
   Panel* panels;  /* Panels stored in this tree */
   unsigned int numpanels;  /* Size of panels */
   Vector3D* points;  /* Points stored in this tree */
   Vector3D* normals;  /* Normals stored in this tree -- added 9/21/08 by JPB for discret. paper*/
   unsigned int numpoints;  /* Size of points, also size of normals */
   unsigned int maxpanelsperfinestcube;  /* How many panels can be in a cube on the finest level */
   unsigned int partitioningdepth;  /* How many levels of cubes */
   BEMKernelType kerneltype;  /* The kernel */
   void* parameters;  /* Parameters passed to Green's function and integrator */
   unsigned int gridpoints;  /* gridpoints^3 is the minimum FFT grid size */
   unsigned int* gridpointsperlevel;  /* number of gridpoints actually used per level >= gridpoints */
   real (*quadraturepoints)[4];  /* Quadrature points for equivalent densities */
   unsigned int numquadraturepoints;  /* Number of quadrature points */
#ifdef POLYNOMIAL
   Matrix* pinv_Fgridpoints;  /* Pseudoinverse of grid points to polynomial coefficients, stored once per level */
#else
   Matrix* pinv_Q2Sgridpoints;  /* Pseudoinverse of grid points to equivalent density sphere quadrature points, stored once per level */
#endif
   real epsilon;  /* SVD precision, relative tolerance */
   BEMLayerType layertype;  /* Single or double layer is needed */
   ComplexSVector**** Tprecomputed;  /* Precomputed FFTs for this tree, indexed by level and x,y,z offsets */
   Cube root;  /* Root node of the cube hierarchy */
   real minimumcubesize;  /* Smallest allowable dimension for a cube */
#ifdef SERIALIZE
   FILE* Dfiles_single[64];
   FILE* Dfiles_double[64];
   FILE* VTfiles[64];
   FILE* Ufiles[64];
   FILE* PVfiles_single[64];
   FILE* PVfiles_double[64];
   FILE* UTIfiles[64];
#endif
} _Tree;

typedef _Tree* Tree;

/* Constructors and Destructors */

Tree Tree_allocate(Panel* panels, unsigned int numpanels, Vector3D* points, unsigned int numpoints, unsigned int maxpanelsperfinestcube, BEMKernelType kerneltype, void* parameters, unsigned int gridpoints, real epsilon, BEMLayerType layertype, real minimumcubesize);
void Tree_free(Tree tree);

/* Operations */

void Tree_lists(Tree tree);
void Tree_fill(Tree tree);
void Tree_memory(Tree tree);
void Tree_multiply(Vector b, Tree tree, Vector x, BEMLayerType layertype);
void Tree_multiply_transpose(Vector b, Tree tree, Vector x, BEMLayerType layertype);
void Tree_multiplyboth(Vector b, Tree tree, Vector xs, Vector xd);
void Tree_writematlabfile(char* filename, Tree tree, BEMLayerType layertype);
void Tree_extractdiagonal(Vector d, Tree tree, BEMLayerType layertype);
