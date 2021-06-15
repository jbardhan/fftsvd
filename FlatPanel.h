/* Typedef */

typedef struct _FlatPanel {
   Vector3D vertex[3];  /* Vertex coordinates in global frame */
   Vector3D centroid;  /* Centroid coordinates in global frame */
   real edgelength[3];  /* Edge lengths */
   Vector3D panelaxis[3];   /* Panel reference frame */
   Vector3D panelvertex[3];  /* Vertex coordinates in panel reference frame */
   real contributionC[3], contributionS[3];  /* cos/sin contribution terms */

   // until we figure out new formulas for contribution[CS]
   Vector3D panelaxisnum[3];
   Vector3D panelvertexnum[3];
   Vector3D edges[3];   /* equations for edges in panel ref frame JPB 11/20/04*/
   real edgeRHS[3]; /* added by JPB for more stable panel transforms */
   real edgeOV[3]; /* vector added by JPB--tells us whether opp. vert is inside*/

   real area;  /* Panel area */
   real moments[16];
   real max_diag;
   real min_diag;

   unsigned int numdirectquadpoints;
   Vector3D *directquadpoints;
   Vector3D *directquadnormals;
   real *directquadweights;
} _FlatPanel;

typedef _FlatPanel* FlatPanel;

/* Constructors and Destructors */

FlatPanel FlatPanel_allocate(Vector3D v1, Vector3D v2, Vector3D v3);
void FlatPanel_free(FlatPanel panel);

void FlatPanel_moments(FlatPanel panel);

/* Operations */

unsigned int FlatPanel_memory(unsigned int order, unsigned int doublelayer);
void FlatPanel_getquadrature(FlatPanel panel, unsigned int rule, unsigned int* numquadpoints,
									  Vector3D** quadpoints, Vector3D** quadnormals, Vector* quadweights);
