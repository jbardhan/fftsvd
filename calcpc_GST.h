#ifndef __CALCPC_GST_H__
#define __CALCPC_GST_H__
#include "GST.h"
void findIntersection(Vector3D point, Vector3D v1, Vector3D v2, real theta);
void Integration_general_GST_single(Vector3D point, GST gst, BEMKernelType kernel, void* parameters, real* slp);
void Integration_general_GST_double(Vector3D point, GST gst, BEMKernelType kernel, void* parameters, real* dlp);
void Integration_oneoverr_GST(Vector3D point, GST gst, void* parameters, real* slp, real *dlp);
void Integration_oneoverr_GST_single(Vector3D point, GST gst, void* parameters, real* slp);
void Integration_oneoverr_GST_double(Vector3D point, GST gst, void* parameters, real* dlp);
void greenInt_XinWang(real *Integrals, FlatRefPanel flatpanel, Vector3D localpnt, unsigned int order);
void Integration_general_GST_single(Vector3D point, GST gst, BEMKernelType kernel, void* parameters, real *slp);
void GST_computePolynomialIntegrals();
void GST_computeVandermondeMatrix();
void GST_computePolynomialValues(); 
real findJacobianDet(Vector3D center, real radius, Vector3D pnt);
void ellipse_int2(real *ellipse_integrals, GST gst, Vector3D point, unsigned int numcoeffs,
						unsigned int* kpvec, unsigned int* epvec,
						BEMKernelType kernel, BEMLayerType layertype, void* parameters);
void conic_section_integrate(real *ellipse_integrals, unsigned int arcnumber,
									  GST gst, Vector3D point, unsigned int numcoeffs,
									  unsigned int* kpvec, unsigned int* epvec,
									  	BEMKernelType kernel, BEMLayerType layertype, void *parameters);
#endif
