Matrix EquivDensity_Q2S_points(Vector3D center, real boxlength, Vector3D* points, unsigned int numpoints, real (*quadraturepoints)[4], unsigned int numquadraturepoints, BEMKernelType kerneltype, void* parameters);
Matrix EquivDensity_Q2S_points_deriv(Vector3D center, real boxlength, Vector3D* points, Vector3D* normals, unsigned int numpoints, real (*quadraturepoints)[4], unsigned int numquadraturepoints, BEMKernelType kerneltype, void* parameters);
Matrix EquivDensity_Q2S_panels(Cube cube, real boxlength, Panel* panels, real (*quadraturepoints)[4], unsigned int numquadraturepoints, BEMKernelType kerneltype, void* parameters, BEMLayerType layertype);

extern real quadrature25[25][4];
extern real quadrature49[49][4];
extern real quadrature64[64][4];
extern real quadrature81[81][4];
extern real quadrature100[100][4];
extern real quadrature121[121][4];

