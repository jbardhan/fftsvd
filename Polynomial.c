#include "FFTSVD.h"

#ifdef POLYNOMIAL

Matrix Polynomial_F_points(Vector3D center, real boxlength, Vector3D* points, unsigned int numpoints, unsigned int gridpoints) {
   unsigned int numcoeffs = gridpoints*(gridpoints+1)*(gridpoints+2)/6;
   unsigned int p, m, i, j;
   real gridspacing = boxlength / (gridpoints - 1);
   Matrix F = Matrix_allocate(numcoeffs, numpoints);

   for (p = 0; p < numpoints; p++) {
      unsigned int count = 0;
      real px = (points[p]->x - center->x) / gridspacing;
      real py = (points[p]->y - center->y) / gridspacing;
      real pz = (points[p]->z - center->z) / gridspacing;     

      for (m = 0; m < gridpoints; m++) {
         for (i = 0; i <= m; i++) {
            for (j = 0; j <= (m-i); j++) {
               F[count][p] = intpow(px, i) * intpow(py, j) * intpow(pz, m-i-j);
               count++;
            }
         }
      }
   }

   return F;
}

Matrix Polynomial_F_panels(Cube cube, real boxlength, Panel* panels, unsigned int gridpoints, BEMLayerType layertype) {
   unsigned int numcoeffs = gridpoints*(gridpoints+1)*(gridpoints+2)/6;
   unsigned int p, m, i, j;
   real gridspacing = boxlength / (gridpoints - 1);
   Matrix F = Matrix_allocate(numcoeffs, cube->numpanelindices);
   Vector3D junkpoint = Vector3D_allocate();
   real parameters[7];

   junkpoint->x = junkpoint->y = junkpoint->z = 1000000.0;

   parameters[3] = cube->center->x;
   parameters[4] = cube->center->y;
   parameters[5] = cube->center->z;
   parameters[6] = gridspacing;

   for (p = 0; p < cube->numpanelindices; p++) {
      unsigned int count = 0;

      for (m = 0; m < gridpoints; m++) {
         for (i = 0; i <= m; i++) {
            for (j = 0; j <= (m-i); j++) {
               parameters[0] = i;
               parameters[1] = j;
               parameters[2] = m-i-j;
               F[count][p] = Integration(junkpoint, panels[cube->panelindices[p]], MONOMIAL_KERNEL, parameters, layertype);
               count++;
            }
         }
      }
   }

   return F;
}

#endif
