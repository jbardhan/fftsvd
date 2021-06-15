real GreensFunction(Vector3D dest, Vector3D src, BEMKernelType kerneltype, void* parameters);
real GreensFunction_deriv(Vector3D dest, Vector3D src, BEMKernelType kerneltype, void* parameters, Vector3D direction);
real GreensFunction_doublederiv(Vector3D dest, Vector3D src, BEMKernelType kerneltype,
										  void* parameters, Vector3D direction, Vector3D direction2);

real GreensFunction_oneoverr(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_oneoverr_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_ekroverr(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_ekroverr_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_ekroverr_desingularized(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_ekroverr_desingularized_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_LJ(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_LJ12(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_LJ6(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_LJ12_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_LJ6_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_monomial(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_monomial_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_inversepower(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_inversepower_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);

// generalized-born kernels
real GreensFunction_Ghosh(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_Ghosh_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_XGB(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_XGB_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_Lesyng(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_Lesyng_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_Grycuk(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_Grycuk_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);

