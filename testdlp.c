#include "FFTSVD.h"

/*#define CONVERSION 332.06*/
#define CONVERSION 332.112

extern unsigned int allcount;
extern unsigned int quadcount;

real maltquad(FlatPanel panel, Panel realpanel) {
   static QuadratureRule qr = NULL;

   if (qr == NULL)
      qr = QuadratureRule_allocate(4);

   unsigned int center, other1, other2;
   unsigned int c, o1, o2;

   Vector3D axis1 = Vector3D_allocate();
   real axis1length;
   Vector3D axis1norm = Vector3D_allocate();
   Vector3D otherside = Vector3D_allocate();
   real axis1dot;

   for (c = 0; c < 3; c++) {
      for (o1 = 0; o1 < 3; o1++) {
         if (o1 == c) continue;

         o2 = 3 - o1 - c;

         Vector3D_sub(axis1, panel->vertex[o1], panel->vertex[c]);
         axis1length = Vector3D_length(axis1);
         Vector3D_copy(axis1norm, axis1);
         Vector3D_scale(axis1norm, 1.0 / axis1length);
         Vector3D_sub(otherside, panel->vertex[o2], panel->vertex[c]);
         axis1dot = Vector3D_dot(otherside, axis1) / Vector3D_dot(axis1, axis1);
         if ((axis1dot > 0.0) && (axis1dot < 1.0)) {
            center = c;
            other1 = o1;
            other2 = o2;
            break;
         }
      }

      if ((axis1dot > 0.0) && (axis1dot < 1.0)) break;
   }

   Vector3D axis1proj = Vector3D_allocate();
   Vector3D_addscaled(axis1proj, panel->vertex[center], axis1dot, axis1);
   Vector3D axis2 = Vector3D_allocate();
   Vector3D_sub(axis2, panel->vertex[other2], axis1proj);
   real height = Vector3D_length(axis2);
   real length1 = Vector3D_distance(axis1proj, panel->vertex[center]);
   real length2 = Vector3D_distance(panel->vertex[other1], axis1proj);
   real slope1 = -length1 / height;
   real slope2 = length2 / height;

   unsigned int i, j;
   Vector3D currentx = Vector3D_allocate();
   Vector3D currentpoint = Vector3D_allocate();

   real integral = 0.0;
   real weightsum = 0.0;

   for (i = 0; i < qr->order; i++) {
      Vector3D_addscaled(currentx, axis1proj, qr->x[i], axis2);
      real starty = slope2 * qr->x[i] * height - length2;
      real endy = slope1 * qr->x[i] * height + length1;
      real dy = endy - starty;
      for (j = 0; j < qr->order; j++) {
         Vector3D_addscaled(currentpoint, currentx, qr->x[j] * dy - endy, axis1norm);
         real currentweight = height * qr->w[i] * dy * qr->w[j];
         real value = 0.0;

         value = Integration(currentpoint, realpanel, POISSON_KERNEL, NULL, DOUBLE_LAYER_INT);

         integral += value * currentweight;
         weightsum += currentweight;
      }
   }

   if (fabs(weightsum - panel->area) > 1e-6)
      printf("MALTQUAD FUTZ\n");

   Vector3D_free(currentpoint);
   Vector3D_free(currentx);
   Vector3D_free(axis2);
   Vector3D_free(axis1proj);
   Vector3D_free(otherside);
   Vector3D_free(axis1norm);
   Vector3D_free(axis1);

   return integral;
}

real galerkin(Panel dest, Panel src) {
   real integral = 0.0;
   unsigned int i;

   for (i = 0; i < src->numdirectquadpoints; i++)
      integral += Integration(src->directquadpoints[i], dest, POISSON_KERNEL, NULL, DOUBLE_LAYER_INT) *
                  src->directquadweights[i];

/*
   if (dest == src)
   integral = 2.0 * M_PI * dest->area;
   else
   integral = maltquad((FlatPanel)(dest->realpanel), src);
*/
   return integral;
}

real split_qual(Vector3D point, GST gst, unsigned int level) {
   Vector3D v01 = Vector3D_allocate();
   Vector3D v12 = Vector3D_allocate();
   Vector3D v20 = Vector3D_allocate();
   Vector3D d = Vector3D_allocate();
   real radius, curdist, integral;
   unsigned int i;

   Vector3D_add(v01, gst->vertices[0], gst->vertices[1]);
   Vector3D_scale(v01, 0.5);
   radius = Vector3D_distance(gst->vertices[0], gst->arccenters[0]);
   Vector3D_sub(d, v01, gst->arccenters[0]);
   Vector3D_addscaled(v01, v01, radius / Vector3D_length(d), d);

   Vector3D_add(v12, gst->vertices[1], gst->vertices[2]);
   Vector3D_scale(v12, 0.5);
   radius = Vector3D_distance(gst->vertices[1], gst->arccenters[1]);
   Vector3D_sub(d, v12, gst->arccenters[1]);
   Vector3D_addscaled(v12, v12, radius / Vector3D_length(d), d);

   Vector3D_add(v20, gst->vertices[2], gst->vertices[0]);
   Vector3D_scale(v20, 0.5);
   radius = Vector3D_distance(gst->vertices[2], gst->arccenters[2]);
   Vector3D_sub(d, v20, gst->arccenters[2]);
   Vector3D_addscaled(v20, v20, radius / Vector3D_length(d), d);

   GST g1 = GST_allocate(gst->center, gst->radius, gst->pit, gst->cavity,
                         gst->vertices[0], v01, v20,
                         gst->arccenters[0], gst->center, gst->arccenters[2]);
   GST g2 = GST_allocate(gst->center, gst->radius, gst->pit, gst->cavity,
                         v01, gst->vertices[1], v12,
                         gst->arccenters[0], gst->arccenters[1], gst->center);
   GST g3 = GST_allocate(gst->center, gst->radius, gst->pit, gst->cavity,
                         v20, v12, gst->vertices[2],
                         gst->center, gst->arccenters[1], gst->arccenters[2]);
   GST g4 = GST_allocate(gst->center, gst->radius, gst->pit, gst->cavity,
                         v01, v12, v20,
                         gst->center, gst->center, gst->center);

   real d1 = 0.0, d2 = 0.0, d3 = 0.0, d4 = 0.0;
   real a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a = 0.0;

   Integration_oneoverr_GST_double(point, g1, NULL, &d1);
   Integration_oneoverr_GST_double(point, g2, NULL, &d2);
   Integration_oneoverr_GST_double(point, g3, NULL, &d3);
   Integration_oneoverr_GST_double(point, g4, NULL, &d4);
/*
   for (i = 0; i < g1->numdirectquadpoints; i++) {
      a1 += g1->directquadweights[i];
      a2 += g2->directquadweights[i];
      a3 += g3->directquadweights[i];
      a4 += g4->directquadweights[i];
      a += gst->directquadweights[i];
   }

   printf("ACHECK: %e %e %e\n", a, a1+a2+a3+a4, fabs(a-a1-a2-a3-a4)/a);
*/
   if (level == 0)
      integral = d1+d2+d3+d4;
   else
      integral = split_qual(point, g1, level+1) + split_qual(point, g2, level+1) +
                 split_qual(point, g3, level+1) + split_qual(point, g4, level+1);

   GST_free(g1);
   GST_free(g2);
   GST_free(g3);
   GST_free(g4);
   Vector3D_free(v01);
   Vector3D_free(v12);
   Vector3D_free(v20);

   return integral;
}

int main(int argc, char* argv[]) {
   FILE* gstfile = NULL;
   FILE* torfile = NULL;
   FILE* vertfile = NULL;
   FILE* facefile = NULL;
   FILE* chargefile = NULL;
   GST* gstpanels;
   TOR* torpanels;
   VertFace vertface;
   Panel* panels;
   unsigned int numgstpanels = 0, numtorpanels = 0, numpanels = 0;
   unsigned int i, j;
   Charge charge;

   if (argc != 4) {
      printf("Usage: %s [gstfile] [torfile] [charge]\n", argv[0]);
      return -1;
   }
/*
   vertface = VertFace_allocate();

   vertfile = fopen(argv[1], "r");

   if (vertfile == NULL) {
      perror("Error opening vertices file");
      return -2;
   }

   VertFace_readvert(vertface, vertfile);

   fclose(vertfile);

   facefile = fopen(argv[2], "r");
    
   if (facefile == NULL) {
      perror("Error opening face file");
      return -2;
   }
    
   VertFace_readface_flip(vertface, facefile);
    
   fclose(facefile);

   VertFace_fix(vertface, 0);

   numpanels = vertface->numfaces;

   //printf("NUMPANELS: %u\n", numpanels);

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   VertFace_getpanels(vertface, panels);
*/

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

   charge = Charge_allocate();

   chargefile = fopen(argv[3], "r");

   if (chargefile == NULL) {
      perror("Error opening charge file");
      return -2; 
   }

   Charge_read(charge, chargefile);

   numpanels = numgstpanels + numtorpanels;

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   for (i = 0; i < numgstpanels; i++) {
      panels[i] = Panel_allocate();
      Panel_GST(panels[i], gstpanels[i], 0);
   }

   for (i = 0; i < numtorpanels; i++) {
      panels[i+numgstpanels] = Panel_allocate();
      Panel_TOR(panels[i+numgstpanels], torpanels[i], 0);
   }

   for (i = 0; i < charge->numcharges; i++) {
      for (j = 0; j < numpanels; j++)
         printf("%g ", 1.0 / Vector3D_distance(charge->points[i], panels[j]->centroid));

      printf("\n");
   }

/*
   FILE* q = fopen("qual", "w");

   for (i = 0; i < numpanels; i++) {
      for (j = 0; j < numpanels; j++) {
        if (i != j) continue;

         real i1 = Integration(panels[j]->centroid, panels[i], POISSON_KERNEL, NULL, DOUBLE_LAYER_INT);
         i1 *= panels[j]->area;

         real i2 = split_qual(panels[j]->centroid, (GST)(panels[i]->realpanel), 0);
         i2 *= panels[j]->area;

         printf("SPLIT: %e %e %e\n", i1, i2, fabs(i1-i2)/fabs(i2));

         fprintf(q, "%e ", i1);
      }

      fprintf(q, "\n");
   }

   fclose(q);
return 0;
   FILE* g = fopen("galerkin", "w");

   for (i = 0; i < numpanels; i++) {
      for (j = 0; j < numpanels; j++) {
         real i2 = galerkin(panels[j], panels[i]);

         fprintf(g, "%e ", i2);
      }

      fprintf(g, "\n");
   }
*/

   return 0;
}
