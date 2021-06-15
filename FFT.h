void FFT_calcDiagonalTranslationOperator(unsigned int level, unsigned int n, Vector3D shift, real gridspacing, BEMKernelType kerneltype, void* parameters, ComplexSVector translate);
void FFT_forwardGridTransform(unsigned int level, unsigned int n, SVector gridCharges, ComplexSVector paddedFourierCharges);
void FFT_backwardGridTransform(unsigned int level, unsigned int n, ComplexSVector paddedFourierRespose, SVector gridPotentials);
