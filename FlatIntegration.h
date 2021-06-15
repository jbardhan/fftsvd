void FlatIntegration_oneoverr_grad(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_oneoverr_deriv(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_oneoverr_deriv_numerical_qual(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_oneoverr_deriv_qual(FlatPanel point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_oneoverr_deriv_qual_point(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);

void FlatIntegration_oneoverr(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_oneoverr_numerical(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_ekroverr_numerical(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_ekroverr_deriv_numerical(Vector3D point, FlatPanel panel, Vector3D normal, void* parameters, real* slp, real* dlp);
void FlatIntegration_ekroverr_desingularized(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_LJ(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_LJ12(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_LJ6(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_Ghosh(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_Lesyng(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_Grycuk(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_general(Vector3D point, FlatPanel panel, BEMKernelType kernel, BEMLayerType layer, void* parameters, real* integral);

real FlatIntegration_maltquad(Vector3D point, FlatPanel panel, BEMKernelType kerneltype, void* parameters, BEMLayerType layertype);

void gen_Stroud_Rule(real *xtab, real *ytab, real *wtab);
