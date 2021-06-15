#ifndef __CALCPC_TOR_H__
#define __CALCPC_TOR_H__
#include "TOR.h"
real realabs(real a);
real callGreensFunction(Vector3D obspnt, Vector3D srcpnt, Vector3D normal,
								BEMKernelType kernel, void *parameters);
void Integration_general_TOR_single(Vector3D point, TOR tor, BEMKernelType kernel, void* parameters, real* slp);
void Integration_oneoverr_TOR(Vector3D point, TOR tor, void* parameters, real* slp, real* dlp);
void Integration_oneoverr_TOR_single(Vector3D point, TOR tor, void* parameters, real* slp);
void Integration_oneoverr_TOR_double(Vector3D point, TOR tor, void* parameters, real* dlp);
void Integration_oneoverr_TOR_double_DJW(Vector3D point, TOR tor, void* parameters, real* dlp);
void Integration_oneoverr_TOR_double_local(Vector3D point, TOR tor, void* parameters, real* dlp);
void Integration_general_TOR_single_local(Vector3D point, TOR tor,
														BEMKernelType kernel, void* parameters, real* slp);
void Integration_general_TOR_single_self(Vector3D point, TOR tor,
													  BEMKernelType kernel, void* parameters, real* slp);
void Integration_general_TOR_single_near(Vector3D point, TOR tor,
													  BEMKernelType kernel, void* parameters, real* slp);

extern unsigned int scount;
extern unsigned int ncount;
extern unsigned int fcount;
#endif
