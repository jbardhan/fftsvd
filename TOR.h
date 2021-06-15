#include "FFTSVD.h"
#ifndef __TOR_H__
#define __TOR_H__

typedef struct _TOR {
  Vector3D centroid;
  Vector3D normalAtCentroid;
  int cavity;
  real area;
  real maxEdgelength;
  real c, a;
  real startPsi, endPsi;
  real startTheta, endTheta;
  Vector3D center; // torus center
  Vector3D normal; // normal of donut plane
  Vector3D localX, localY;
  Vector3D Tmat[3];
  Vector3D Tvec;
  int isLocal; // should be deprecated soon...
  unsigned int thetaIndex, psiIndex;
  unsigned int numthetadivs, numpsidivs;
  unsigned int torusID;

  unsigned int numdirectquadpoints;
  Vector3D *directquadpoints;
  Vector3D *directquadnormals;
  real *directquadweights;
} _TOR;

typedef struct _TOR* TOR;

TOR TOR_allocate(real torusRadius, real probeRadius, real startTheta, real endTheta,
					  real startPsi, real endPsi, Vector3D center, Vector3D normal,
					  Vector3D localX, Vector3D localY, unsigned int thetaIndex, unsigned int psiIndex, int cavity,
					  unsigned int numthetadivs, unsigned int numpsidivs, unsigned int torusID);
void TOR_free(TOR t);
void TOR_readfile(unsigned int *numTORpanels, TOR** torList, FILE* file, unsigned int ignorecav);
void torToCart(Vector3D point, Vector3D normal,
					real c, real a, real theta, real psi);
#endif
