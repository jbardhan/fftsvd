/* Typedef */

typedef enum { FLAT_PANEL, GST_PANEL, TOR_PANEL, POINT_PANEL } PanelType;

typedef struct _Panel {
   Vector3D centroid;  /* Centroid coordinates in global frame */
   Vector3D normal;  /* Normal vector at centroid */
   real area;  /* Panel area */
   real maxedgelength;

   unsigned int numdirectquadpoints;
   Vector3D* directquadpoints;
   Vector3D* directquadnormals;
   Vector directquadweights;

#ifdef ENABLE_GALERKIN
   unsigned int numgalerkinquadpoints;
   Vector3D* galerkinquadpoints;
   Vector3D* galerkinquadnormals;
   Vector galerkinquadweights;
#endif

   PanelType type;
   void* realpanel;
} _Panel;

typedef _Panel* Panel;

/* Constructors and Destructors */

Panel Panel_allocate();
void Panel_free(Panel);

/* Operators */

void Panel_FlatPanel(Panel panel, FlatPanel flatpanel);
void Panel_GST(Panel Panel, GST gst, unsigned int usecom);
void Panel_TOR(Panel panel, TOR tor, unsigned int usecom);
void Panel_Vector3D(Panel panel, Vector3D v3d);

real Panel_memory(Panel panel, BEMLayerType layertype);

real Panel_quadrature(Vector3D point, Panel panel, BEMKernelType kerneltype, void* parameters, BEMLayerType layertype);
