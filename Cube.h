/* Typedef */

typedef struct _Cube {
   unsigned int level;  /* depth of this cube in the tree */
   unsigned int indices[3];  /* grid coordinates of this cube on this level */
   Vector3D bounds[2];  /* diagonal corners of the cube */
   Vector3D center; /* center of the cube */
   struct _Cube* parent;  /* Parent cube */
   unsigned int leaf;  /* 1 if a leaf, otherwise 0 */
   unsigned int* panelindices;  /* indices of the panels contained in this cube */
   unsigned int numpanelindices;  /* number of panels in this cube */
   unsigned int* pointindices;  /* indices of the collocation points contained in this cube */
   unsigned int numpointindices;  /* number of collocation points in this cube */
   struct _Cube** localcubes;  /* pointers to local cubes */
   unsigned int numlocalcubes;   /* number of local cubes */
   struct _Cube** interactingcubes;  /* pointers to interacting cubes */
   unsigned int numinteractingcubes;  /* number of interacting cubes */
   struct _Cube* children[2][2][2];   /* pointers to children */
   unsigned int drows, dcolumns;   /* dimensions for D if a leaf (numpanelindices + numpanelindices in local cubes) */
   unsigned int rowrank, columnrank;   /* inner dimensions for VTsrc and Udest */
   SMatrix VTsrc;  /* row basis for source panels | interaction points (rowrank x panelindices) for single and double layer */
   SMatrix Udest;  /* column basis for source points | interaction panels (pointindices x columnrank) for single and double layer */
   SMatrix D_single;  /* dense matrix for self+local if a leaf for single layer */
   SMatrix D_double;  /* dense matrix for self+local if a leaf for double layer */
#ifdef ADAPTIVE
   SVector VT_q;  /* rowrank temporary storage for VT*q during matrix-vector multiply */
   SVector sum_K_VT_q;  /* columnrank temporary storage for K*VT*q during matrix-vector multiply */
#endif
   ComplexSVector FFT_PV_VT_q;  /* ((2*gridpoints-1)^3) temporary storage for FFT*PV*VT*q during matrix-vector multiply */
   ComplexSVector sum_T_FFT_PV_VT_q;  /* ((2*gridpoints-1)^3)) temporary storage for sum (T.*FFT*PV*VT*q) during matrix-vector multiply */
   ComplexSVector FFT_UTIT_UT_q;  /* ((2*gridpoints-1)^3) temporary storage for FFT*UTIT*UT*q during matrix-vector transpose multiply */
   ComplexSVector sum_T_FFT_UTIT_UT_q;  /* ((2*gridpoints-1)^3)) temporary storage for sum (T.*FFT*UTIT*UT*q) during matrix-vector transpose multiply */
   SMatrix UTI;  /* UTdest * interpolation (columnrank * gridpoints^3) */
   SMatrix PV_single;  /* projection * Vsrc (gridpoints^3 * rowrank) for single layer */
   SMatrix PV_double;  /* projection * Vsrc (gridpoints^3 * rowrank) for double layer */
   SMatrix P_single;  /* projection (gridpoints^3 * numpanels) for single layer */
   SMatrix P_double;  /* projection (gridpoints^3 * numpanels) for double layer */
   SMatrix I_; /* interpolation (numpoints * gridpoints^3) */
   SMatrix I_double; /* interpolation for potential gradient dotted with normals... (numpoints * gridpoints^3).  from galerkin branch */
#ifdef ADAPTIVE
   SMatrix* K_single;  /* K-matrices for single layer */
   SMatrix* K_double;  /* K-matrices for double layer */
#endif
   struct _Tree* tree; /* the tree this cube belongs to */
#ifdef SERIALIZE
   unsigned int serialid;
#endif
#ifdef SCATTER
  unsigned int isInsideProtein;
  real densityInsideCube;
#endif
} _Cube;

typedef _Cube* Cube;

/* Constructors and Destructors */

Cube Cube_allocate(Panel* panels, unsigned int* panelindices, unsigned int numpanelindices, Vector3D* points, unsigned int* pointindices, unsigned int numpointindices, unsigned int level, unsigned int indices[3], Vector3D bounds[2], Cube parent, struct _Tree* tree, unsigned int partitioningdepth);
void Cube_free(Cube cube);

/* Operations */

void Cube_lists(Cube cube);
void Cube_fill_D(Cube cube);
void Cube_fill_P_I(Cube cube, unsigned int recurse);
void Cube_fill_U_VT(Cube cube);
void Cube_fill_PV_UTI(Cube cube);
#ifdef ADAPTIVE
void Cube_fill_K(Cube cube);
#endif
void Cube_memory(Cube cube, unsigned int* Cubemem, unsigned int* Umem, unsigned int* VTmem, unsigned int* UTImem, unsigned int* PVmem, unsigned int* Kmem, unsigned int* Dmem);
void Cube_multiply_D(Vector b, Cube cube, Vector x, BEMLayerType layertype);
void Cube_multiply_DT(Vector b, Cube cube, Vector x, BEMLayerType layertype);
void Cube_multiply_FFT_PV_VT(Cube cube, Vector x, BEMLayerType layertype);
void Cube_multiply_FFT_UTIT_UT(Cube cube, Vector x);
void Cube_multiplyboth_FFT_PV_VT(Cube cube, Vector xs, Vector xd);
void Cube_multiply_T(Cube cube);
void Cube_multiply_T_transpose(Cube cube);
void Cube_multiply_U_UTI_IFFT(Vector b, Cube cube);
void Cube_multiply_V_PVT_IFFT(Vector b, Cube cube, BEMLayerType layertype);
void Cube_multiply_K(Cube cube, BEMLayerType layertype);
void Cube_clear_tempvectors(Cube cube);
void Cube_leafpanelstats(Cube cube, unsigned int* numleaves, unsigned int* maxpanels, unsigned int* numleaveswithinlimitpanels, unsigned int panellimit);
void Cube_extractdiagonal(Vector diag, Cube cube, BEMLayerType layertype);

#ifdef SCATTER
extern unsigned int trackCube;
#endif
