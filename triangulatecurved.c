#include "FFTSVD.h"

#define MAXLEVELS 6
#define MAXVERTICES 2000000
#define MAXFACES 2000000
#define TOL 1e-12

real colors [4][3] = {{1.0, 0.0, 0.0},
                      {0.0, 1.0, 0.0},
                      {0.0, 0.0, 1.0},
                      {1.0, 1.0, 0.0}};

real (**unitvertices)[2];
unsigned int *numunitvertices;
unsigned int (**unitfaces)[3];
unsigned int *numunitfaces;

Vector3D* vertices;
Vector3D* normals;
unsigned int (*faces)[3];
unsigned int* panels;
unsigned int* type;

real (*edgepoints)[3];

real maxarea = 1.0;

unsigned int vcount = 0, fcount = 0, epcount = 0;

void unittriangulate() {
   unsigned int i;

   unitvertices = (real(**)[2])calloc(MAXLEVELS, sizeof(real(*)[2]));
   numunitvertices = (unsigned int*)calloc(MAXLEVELS, sizeof(unsigned int));
   unitfaces = (unsigned int(**)[3])calloc(MAXLEVELS, sizeof(unsigned int(*)[3]));
   numunitfaces = (unsigned int*)calloc(MAXLEVELS, sizeof(unsigned int));

   numunitvertices[0] = 3;
   unitvertices[0] = (real(*)[2])calloc(numunitvertices[0], sizeof(real[2]));
   numunitfaces[0] = 1;
   unitfaces[0] = (unsigned int(*)[3])calloc(numunitfaces[0], sizeof(unsigned int[3]));

   unitvertices[0][0][0] = 0.0;
   unitvertices[0][0][1] = 0.0;
   unitvertices[0][1][0] = 0.0;
   unitvertices[0][1][1] = 1.0;
   unitvertices[0][2][0] = 1.0;
   unitvertices[0][2][1] = 0.0;

   unitfaces[0][0][0] = 0;
   unitfaces[0][0][1] = 1;
   unitfaces[0][0][2] = 2;
   
   for (i = 1; i < MAXLEVELS; i++) {
      unsigned int v, f, vc, fc;

      numunitvertices[i] = 4 * numunitvertices[i-1] - 6;
      unitvertices[i] = (real(*)[2])calloc(numunitvertices[i], sizeof(real[2]));
      numunitfaces[i] = 4 * numunitfaces[i-1];
      unitfaces[i] = (unsigned int(*)[3])calloc(numunitfaces[i], sizeof(unsigned int[3]));

      for (v = 0; v < numunitvertices[i-1]; v++) {
         unitvertices[i][v][0] = unitvertices[i-1][v][0];
         unitvertices[i][v][1] = unitvertices[i-1][v][1];
      }

      vc = numunitvertices[i-1];
      fc = 0;

      for (f = 0; f < numunitfaces[i-1]; f++) {
         real v1[2], v2[2], v3[2], m12[2], m23[2], m31[2];
         v1[0] = unitvertices[i-1][unitfaces[i-1][f][0]][0];
         v1[1] = unitvertices[i-1][unitfaces[i-1][f][0]][1];
         v2[0] = unitvertices[i-1][unitfaces[i-1][f][1]][0];
         v2[1] = unitvertices[i-1][unitfaces[i-1][f][1]][1];
         v3[0] = unitvertices[i-1][unitfaces[i-1][f][2]][0];
         v3[1] = unitvertices[i-1][unitfaces[i-1][f][2]][1];

         m12[0] = 0.5 * (v1[0] + v2[0]);
         m12[1] = 0.5 * (v1[1] + v2[1]);
         m23[0] = 0.5 * (v2[0] + v3[0]);
         m23[1] = 0.5 * (v2[1] + v3[1]);
         m31[0] = 0.5 * (v3[0] + v1[0]);
         m31[1] = 0.5 * (v3[1] + v1[1]);

         unitvertices[i][vc][0] = m12[0];
         unitvertices[i][vc][1] = m12[1];
         unitvertices[i][vc+1][0] = m23[0];
         unitvertices[i][vc+1][1] = m23[1];
         unitvertices[i][vc+2][0] = m31[0];
         unitvertices[i][vc+2][1] = m31[1];

         unitfaces[i][fc][0] = unitfaces[i-1][f][0];
         unitfaces[i][fc][1] = vc;
         unitfaces[i][fc][2] = vc+2;
         unitfaces[i][fc+1][0] = vc;
         unitfaces[i][fc+1][1] = unitfaces[i-1][f][1];
         unitfaces[i][fc+1][2] = vc+1;
         unitfaces[i][fc+2][0] = vc+1;
         unitfaces[i][fc+2][1] = unitfaces[i-1][f][2];
         unitfaces[i][fc+2][2] = vc+2;
         unitfaces[i][fc+3][0] = vc;
         unitfaces[i][fc+3][1] = vc+1;
         unitfaces[i][fc+3][2] = vc+2;

         vc += 3;
         fc += 4;
      }
   }
}

unsigned int ongstedge(real eta1, real ksi1, real eta2, real ksi2) {
   if ((fabs(eta1) < TOL) && (fabs(eta2) < TOL))
      return 1;
   else if ((fabs(ksi1) < TOL) && (fabs(ksi2) < TOL))
      return 1;
   else if ((fabs(eta1+ksi1-1.0) < TOL) && (fabs(eta2+ksi2-1.0) < TOL))
      return 1;
   else
      return 0;
}

unsigned int ontoredge(real eta1, real ksi1, real eta2, real ksi2) {
   if ((fabs(eta1) < TOL) && (fabs(eta2) < TOL))
      return 1;
   else if ((fabs(ksi1) < TOL) && (fabs(ksi2) < TOL))
      return 1;
   else if ((fabs(eta1 - 1.0) < TOL) && (fabs(eta2 - 1.0) < TOL))
      return 1;
   else if ((fabs(ksi1 - 1.0) < TOL) && (fabs(ksi2 - 1.0) < TOL))
      return 1;
   else
      return 0;
}

void triangulategst(GST gst, unsigned int panel) {
   unsigned int i, next, f, level = 0;
   real area = 0.0, ratio, junk;

   for (i = 0; i < gst->numdirectquadpoints; i++)
      area += gst->directquadweights[i];

   ratio = area / maxarea;

   while (intpow(4, level) < ratio) level++;

   if (level >= MAXLEVELS) {
      printf("Increase MAXLEVELS!\n");
      exit(-3);
   }

   Vector3D dTmat[3];
   Vector3D dTvec;
   Curve edges[3];
   dTmat[0] = Vector3D_allocate();
   dTmat[1] = Vector3D_allocate();
   dTmat[2] = Vector3D_allocate();
   dTvec = Vector3D_allocate();

   GST_getStandardPosition(gst, edges, dTmat, dTvec);

   for (f = 0; f < numunitfaces[level]; f++) {
      faces[fcount][0] = vcount;
      faces[fcount][1] = vcount+1;
      faces[fcount][2] = vcount+2;

      panels[fcount] = panel;
      type[fcount] = gst->pit;

      for (i = 0; i < 3; i++) {
         real eta = unitvertices[level][unitfaces[level][f][i]][0];
         real ksi = unitvertices[level][unitfaces[level][f][i]][1];
         vertices[vcount] = Vector3D_allocate();
         normals[vcount] = Vector3D_allocate();

         if (ksi == 0.0) ksi = TOL;
         if (ksi == 1.0) ksi = 1-TOL;
         if (eta == 0.0) eta = TOL;
         if (eta == 1.0) eta = 1-TOL;

         GST_get_direct_quadPoint(vertices[vcount], normals[vcount], &junk,
                                  eta, ksi, gst, edges, dTmat, dTvec);

         if (gst->pit > 0) 
            Vector3D_scale(normals[vcount], -1.0);
         if (gst->cavity > 0)
            Vector3D_scale(normals[vcount], -1.0);

         if (isnan(vertices[vcount]->x)) {
            printf("Coordinate transform failed\n");
            printf("%g %g\n", eta, ksi);
            exit(-4);
         }

         vcount++;

         if (vcount == MAXVERTICES) {
            printf("Increase MAXVERTICES\n");
            exit(-5);
         }
      }

      for (i = 0; i < 3; i++) {
         next = i + 1;
         if (next == 3) next = 0;

         real eta1 = unitvertices[level][unitfaces[level][f][i]][0];
         real ksi1 = unitvertices[level][unitfaces[level][f][i]][1];
         real eta2 = unitvertices[level][unitfaces[level][f][next]][0];
         real ksi2 = unitvertices[level][unitfaces[level][f][next]][1];

         if (ongstedge(eta1, ksi1, eta2, ksi2)) {
            if (epcount >= MAXVERTICES) {
               printf("Increase MAXVERTICES\n");
               exit(-5);
            }

            edgepoints[epcount][0] = vertices[vcount-3+i]->x;
            edgepoints[epcount][1] = vertices[vcount-3+i]->y;
            edgepoints[epcount][2] = vertices[vcount-3+i]->z;
            edgepoints[epcount+1][0] = vertices[vcount-3+next]->x;
            edgepoints[epcount+1][1] = vertices[vcount-3+next]->y;
            edgepoints[epcount+1][2] = vertices[vcount-3+next]->z;

            epcount += 2;
         }
      }

      fcount++;

      if (vcount == MAXFACES) {
         printf("Increase MAXFACES\n");
         exit(-5);
      }
   }
}

void triangulatetor(TOR tor, unsigned int panel) {
   unsigned int i, next, f, level = 0;
   real area = 0.0, ratio, junk;

   for (i = 0; i < tor->numdirectquadpoints; i++)
      area += tor->directquadweights[i];

   ratio = (0.5 * area) / maxarea;

   while (intpow(4, level) < ratio) level++;

   if (level >= MAXLEVELS) {
      printf("Increase MAXLEVELS!\n");
      exit(-3);
   }

   for (f = 0; f < numunitfaces[level]; f++) {
      faces[fcount][0] = vcount;
      faces[fcount][1] = vcount+1;
      faces[fcount][2] = vcount+2;

      panels[fcount] = panel;
      type[fcount] = 2;

      for (i = 0; i < 3; i++) {
         real eta = unitvertices[level][unitfaces[level][f][i]][0];
         real ksi = unitvertices[level][unitfaces[level][f][i]][1];
         vertices[vcount] = Vector3D_allocate();
         normals[vcount] = Vector3D_allocate();

         real theta = tor->startTheta + eta * (tor->endTheta - tor->startTheta);
         real psi = tor->startPsi + ksi * (tor->endPsi - tor->startPsi);

         torToCart(vertices[vcount], normals[vcount], tor->c, tor->a, theta, psi);

         Vector3D_sub(vertices[vcount], vertices[vcount], tor->Tvec);
         Vector3D_transformVecs_inverse(vertices[vcount], tor->Tmat[0], tor->Tmat[1], tor->Tmat[2], vertices[vcount]);
         Vector3D_transformVecs_inverse(normals[vcount], tor->Tmat[0], tor->Tmat[1], tor->Tmat[2], normals[vcount]);

         Vector3D_scale(normals[vcount], -1.0);
         if (tor->cavity > 0)
            Vector3D_scale(normals[vcount], -1.0);

         if (isnan(vertices[vcount]->x)) {
            printf("Coordinate transform failed\n");
            printf("%g %g\n", eta, ksi);
            exit(-4);
         }

         vcount++;

         if (vcount == MAXVERTICES) {
            printf("Increase MAXVERTICES\n");
            exit(-5);
         }
      }

      for (i = 0; i < 3; i++) {
         next = i + 1;
         if (next == 3) next = 0;

         real eta1 = unitvertices[level][unitfaces[level][f][i]][0];
         real ksi1 = unitvertices[level][unitfaces[level][f][i]][1];
         real eta2 = unitvertices[level][unitfaces[level][f][next]][0];
         real ksi2 = unitvertices[level][unitfaces[level][f][next]][1];

         if (ontoredge(eta1, ksi1, eta2, ksi2)) {
            if (epcount >= MAXVERTICES) {
               printf("Increase MAXVERTICES\n");
               exit(-5);
            }

            edgepoints[epcount][0] = vertices[vcount-3+i]->x;
            edgepoints[epcount][1] = vertices[vcount-3+i]->y;
            edgepoints[epcount][2] = vertices[vcount-3+i]->z;
            edgepoints[epcount+1][0] = vertices[vcount-3+next]->x;
            edgepoints[epcount+1][1] = vertices[vcount-3+next]->y;
            edgepoints[epcount+1][2] = vertices[vcount-3+next]->z;

            epcount += 2;
         }
      }

      fcount++;

      if (vcount == MAXFACES) {
         printf("Increase MAXFACES\n");
         exit(-5);
      }
   }

   for (f = 0; f < numunitfaces[level]; f++) {
      faces[fcount][0] = vcount;
      faces[fcount][1] = vcount+1;
      faces[fcount][2] = vcount+2;

      panels[fcount] = panel;
      type[fcount] = 2;

      for (i = 0; i < 3; i++) {
         real eta = unitvertices[level][unitfaces[level][f][i]][0];
         real ksi = unitvertices[level][unitfaces[level][f][i]][1];
         vertices[vcount] = Vector3D_allocate();
         normals[vcount] = Vector3D_allocate();

         eta = 1.0 - eta;
         ksi = 1.0 - ksi;

         real theta = tor->startTheta + eta * (tor->endTheta - tor->startTheta);
         real psi = tor->startPsi + ksi * (tor->endPsi - tor->startPsi);

         torToCart(vertices[vcount], normals[vcount], tor->c, tor->a, theta, psi);

         Vector3D_sub(vertices[vcount], vertices[vcount], tor->Tvec);
         Vector3D_transformVecs_inverse(vertices[vcount], tor->Tmat[0], tor->Tmat[1], tor->Tmat[2], vertices[vcount]);
         Vector3D_transformVecs_inverse(normals[vcount], tor->Tmat[0], tor->Tmat[1], tor->Tmat[2], normals[vcount]);

         Vector3D_scale(normals[vcount], -1.0);
         if (tor->cavity > 0)
            Vector3D_scale(normals[vcount], -1.0);

         if (isnan(vertices[vcount]->x)) {
            printf("Coordinate transform failed\n");
            printf("%g %g\n", eta, ksi);
            exit(-4);
         }

         vcount++;

         if (vcount == MAXVERTICES) {
            printf("Increase MAXVERTICES\n");
            exit(-5);
         }
      }

      for (i = 0; i < 3; i++) {
         next = i + 1;
         if (next == 3) next = 0;

         real eta1 = unitvertices[level][unitfaces[level][f][i]][0];
         real ksi1 = unitvertices[level][unitfaces[level][f][i]][1];
         real eta2 = unitvertices[level][unitfaces[level][f][next]][0];
         real ksi2 = unitvertices[level][unitfaces[level][f][next]][1];

         if (ontoredge(eta1, ksi1, eta2, ksi2)) {
            if (epcount >= MAXVERTICES) {
               printf("Increase MAXVERTICES\n");
               exit(-5);
            }

            edgepoints[epcount][0] = vertices[vcount-3+i]->x;
            edgepoints[epcount][1] = vertices[vcount-3+i]->y;
            edgepoints[epcount][2] = vertices[vcount-3+i]->z;
            edgepoints[epcount+1][0] = vertices[vcount-3+next]->x;
            edgepoints[epcount+1][1] = vertices[vcount-3+next]->y;
            edgepoints[epcount+1][2] = vertices[vcount-3+next]->z;

            epcount += 2;
         }
      }

      fcount++;

      if (vcount == MAXFACES) {
         printf("Increase MAXFACES\n");
         exit(-5);
      }
   }
}

void writeoff(FILE* offfile) {
   unsigned int i;
   real r, g, b;
   int currentpanel = -1;

   fprintf(offfile, "NOFF\n");
   fprintf(offfile, "%u %u 0\n", vcount, fcount);

   for (i = 0; i < vcount; i++)
      fprintf(offfile, "%f %f %f %f %f %f\n", vertices[i]->x, vertices[i]->y, vertices[i]->z, normals[i]->x, normals[i]->y, normals[i]->z);
   
   for (i = 0; i < fcount; i++) {
/*
      r = colors[i % 4][0];
      g = colors[i % 4][1];
      b = colors[i % 4][2];
*/

      //r = 1.0; b = 0.0; g = 0.0;
      r = colors[type[i]][0];
      g = colors[type[i]][1];
      b = colors[type[i]][2];
/*
      if (currentpanel != panels[i]) {
         r = (real)random() / (real)RAND_MAX;
         g = (real)random() / (real)RAND_MAX;
         b = (real)random() / (real)RAND_MAX;
         currentpanel = panels[i];
      }
*/
      fprintf(offfile, "3 %u %u %u %f %f %f\n", faces[i][0], faces[i][1], faces[i][2], r, g, b);
   }
}

void writevect(FILE* vectfile) {
   unsigned int i;

   fprintf(vectfile, "VECT\n");
   fprintf(vectfile, "%u %u 0\n", epcount/2, epcount);
   for (i = 0; i < epcount/2; i++)
      fprintf(vectfile, "2 ");
   fprintf(vectfile, "\n");
   for (i = 0; i < epcount/2; i++)
      fprintf(vectfile, "0 ");
   fprintf(vectfile, "\n");

   for (i = 0; i < epcount; i+=2)
      fprintf(vectfile, "%f %f %f %f %f %f\n", edgepoints[i][0], edgepoints[i][1], edgepoints[i][2], edgepoints[i+1][0], edgepoints[i+1][1], edgepoints[i+1][2]);
}

int main(int argc, char* argv[]) {
   FILE* gstfile = NULL;
   FILE* torfile = NULL;
   FILE* offfile = NULL;
   FILE* vectfile = NULL;
   GST* gstpanels;
   TOR* torpanels;
   unsigned int numgstpanels = 0, numtorpanels = 0, i;

   if (argc != 6) {
      printf("Usage: %s [gstfile] [torfile] [offfile] [vectfile] [maxarea]\n", argv[0]);
      return -1;
   }

   gstfile = fopen(argv[1], "r");

   if (gstfile == NULL) {
      perror("Error opening GST file");
      return -2;
   }

   GST_readfile(&numgstpanels, &gstpanels, gstfile, 0);

   fclose(gstfile);

   torfile = fopen(argv[2], "r");
    
   if (torfile == NULL) {
      perror("Error opening TOR file");
      return -2;
   }

   TOR_readfile(&numtorpanels, &torpanels, torfile, 0);
    
   fclose(torfile);

   maxarea = atof(argv[5]);

   unittriangulate();

   vertices = (Vector3D*)calloc(MAXVERTICES, sizeof(Vector3D));
   normals = (Vector3D*)calloc(MAXVERTICES, sizeof(Vector3D));
   faces = (unsigned int(*)[3])calloc(MAXFACES, sizeof(unsigned int[3]));
   panels = (unsigned int*)calloc(MAXFACES, sizeof(unsigned int));
   type = (unsigned int*)calloc(MAXFACES, sizeof(unsigned int));
   edgepoints = (real(*)[3])calloc(MAXVERTICES, sizeof(real[3]));

   for (i = 0; i < numgstpanels; i++)
      triangulategst(gstpanels[i], i);

   for (i = 0; i < numtorpanels; i++)
      triangulatetor(torpanels[i], numgstpanels+i);

   offfile = fopen(argv[3], "w");

   if (offfile == NULL) {
      perror("Error opening off file");
      return -2;
   }

   writeoff(offfile);

   fclose(offfile);

   vectfile = fopen(argv[4], "w");

   if (vectfile == NULL) {
      perror("Error opening vect file");
      return -2;
   }

   writevect(vectfile);

   fclose(vectfile);

   return 0;
}
