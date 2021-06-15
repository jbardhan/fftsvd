#include "FFTSVD.h"
#ifndef __GST_H__
#define __GST_H__
typedef struct _FlatRefPanel {
  Vector3D centroid;
  Vector3D normal;
  real rhs;
  Vector3D vertices[3];
  unsigned int numquadpoints;
  real *ksi;
  real *eta;
  real *weights;
  real maxEdgelength;
} _FlatRefPanel;
typedef _FlatRefPanel* FlatRefPanel;

FlatRefPanel FlatRefPanel_allocate(Vector3D v1, Vector3D v2, Vector3D v3,
									  Vector3D centroid, Vector3D normal,
									  real rhs, unsigned int numquadpoints,
									  real *ksi, real *eta, real *weights,
									  real maxEdgelength);
void FlatRefPanel_free(FlatRefPanel fp);

typedef struct _Conic {
  real asquared;
  real bsquared;
  real theta1Int;
  real theta2Int;
  Vector3D Tmat[3];
  Vector3D Tvec;
} _Conic;

typedef _Conic* Conic;

Conic Conic_allocate(real asq, real bsq, real theta1Int, real theta2Int,
							Vector3D *Tmat, Vector3D Tvec);

void Conic_free(Conic c);

typedef struct _GST {
  Vector3D center;
  Vector3D vertices[3];
  Vector3D arccenters[3];
  int pit, cavity;
  real radius;
  real area;

  Vector3D Tmat[3];
  Vector3D Tvec;

  Vector3D centroid; // on GST
  Vector3D normalAtCentroid; 
  FlatRefPanel flatpanel;
  int factor[3];  // add or subtract or zero
  Conic conic[3];

  unsigned int numdirectquadpoints;
  Vector3D *directquadpoints;
  Vector3D *directquadnormals;
  real *directquadweights;
  
  Matrix vandermondeInverse;
} _GST;
typedef _GST* GST;

typedef struct _Curve {
  Vector3D ac;
  Vector3D v1;
  Vector3D v2;
  Vector3D X; // points from ac -> v1
  Vector3D Y;

  real radius;
  real endTheta;
} _Curve;

typedef struct _Curve* Curve;
Curve Curve_allocate(Vector3D ac, Vector3D v1, Vector3D v2);
void Curve_getParamPoint(Vector3D p, Vector3D dpdalpha, Curve c, real alpha);
void Curve_free(Curve c);
void PermuteVectors(Vector3D v0, Vector3D v1, Vector3D v2);

GST GST_allocate(Vector3D center, real radius, int caporpit, int cavity,
					  Vector3D v1, Vector3D v2, Vector3D v3,
					  Vector3D ac1, Vector3D ac2, Vector3D ac3);
void GST_free(GST gst);

void GST_readfile(unsigned int *numGSTpanels, GST** gstList, FILE* file, unsigned int ignorecav);
void findSpherePoint(Vector3D sphvert, Vector3D center,
							real radius, Vector3D oldvert);
void initialize_GST_FlatRefPanel(GST gst);
void generateQuadraturePoints(Vector3D* vertices, unsigned int *numquadpoints, real **ksi, real **eta, real **weights);
void figure_out_correction_sign(GST gst, unsigned int index);
void conic_section_identify(GST gst, unsigned int index);
void projectToPlane(Vector3D onplane, Vector3D start, Vector3D end);
void GST_get_conic(real* Vret, Vector3D *M, unsigned int numpoints);
void GST_get_direct_quadPoint(Vector3D point, Vector3D normal, real *detJ,
										real ksi, real eta, GST gst,
										Curve *edges, Vector3D *Tmat, Vector3D Tvec);

void GST_getStandardPosition(GST gst, Curve *edges, Vector3D *Tmat, Vector3D Tvec);
void GST_sphGSTtoCart(real radius, real theta, real phi, Vector3D point, Vector3D normal, Vector3D dPoint_dTheta, Vector3D dPoint_dPhi);
void GST_getCircleArcIntersection(Vector3D center, Vector3D normal, real radius,
											 Curve curve, Vector3D intPoint,
											 real *alpha, Vector3D dIntPoint_dAlpha);
void GST_getCentroid(GST gst);
#endif
