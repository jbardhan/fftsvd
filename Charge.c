#include "FFTSVD.h"

/* Constructors and Destructors */

Charge Charge_allocate() {
   return (Charge)calloc(1, sizeof(_Charge));
}

void Charge_free(Charge charge) {
   unsigned int i;
   
   for (i = 0; i < charge->numcharges; i++)
      Vector3D_free(charge->points[i]);
   
   free(charge->points);
   Vector_free(charge->charges);

   free(charge);
}

/* Operations */

void Charge_read(Charge charge, FILE* file) {
   char line[128];
   unsigned int i;
        
   charge->numcharges = 0;
    
   while (fgets(line, 128, file))
      charge->numcharges++;
       
   charge->points = (Vector3D*)calloc(charge->numcharges, sizeof(Vector3D));
   charge->charges = Vector_allocate(charge->numcharges);
    
   rewind(file);
    
   for (i = 0; i < charge->numcharges; i++) {
      charge->points[i] = Vector3D_allocate();

      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf %lf", &charge->points[i]->x, &charge->points[i]->y, &charge->points[i]->z, &charge->charges[i]);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f %f", &charge->points[i]->x, &charge->points[i]->y, &charge->points[i]->z, &charge->charges[i]);
#endif
#endif
   }
}

void Charge_read_centroids(Charge charge, FILE* file, Vector3D* centroids, unsigned int numcentroids) {
   char line[128];
   unsigned int i, numlines = 0;;
        
   charge->numcharges = numcentroids;
    
   while (fgets(line, 128, file)) {
      numlines++;
      charge->numcharges++;
   }
       
   charge->points = (Vector3D*)calloc(charge->numcharges, sizeof(Vector3D));
   charge->charges = Vector_allocate(charge->numcharges);
    
   rewind(file);
    
   for (i = 0; i < numlines; i++) {
      charge->points[i] = Vector3D_allocate();

      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf %lf", &charge->points[i]->x, &charge->points[i]->y, &charge->points[i]->z, &charge->charges[i]);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f %f", &charge->points[i]->x, &charge->points[i]->y, &charge->points[i]->z, &charge->charges[i]);
#endif
#endif
   }

   for (i = 0; i < numcentroids; i++) {
      charge->points[i+numlines] = Vector3D_allocate();
      Vector3D_copy(charge->points[i+numlines], centroids[i]);
      charge->charges[i+numlines] = 0.0;
   }
}

void Charge_makerhs(Vector rhs, Charge charge, Panel* panels, unsigned int numpanels) {
   unsigned int i, j;
    
   for (i = 0; i < numpanels; i++) {
      rhs[i] = (real)0;
      for (j = 0; j < charge->numcharges; j++)
         rhs[i] += charge->charges[j] / Vector3D_distance(charge->points[j], panels[i]->centroid);
   }
}

void Charge_makerhs_ecf(Vector rhs, Charge charge, Panel* panels, unsigned int numpanels) {
   unsigned int i, j;
    
   for (i = 0; i < numpanels; i++) {
      rhs[i] = (real)0;
      for (j = 0; j < charge->numcharges; j++) {
         Vector3D rvector = Vector3D_allocate();
         Vector3D_sub(rvector, panels[i]->centroid, charge->points[j]);
         real r = Vector3D_length(rvector);
         Vector3D_normalize(rvector);

         rhs[i] += charge->charges[j] * (80.0 - 4.0) / (4.0 * M_PI * 4.0) / (r*r) * Vector3D_dot(rvector, panels[i]->normal);

         Vector3D_free(rvector);
      }
   }
}

void Charge_makerhs_juffersimple(Vector rhs, Charge charge, Panel* panels, unsigned int numpanels, real idiel, real odiel) {
  unsigned int i, j;
  real r;
  Vector3D rvector = Vector3D_allocate();
  for (i = 0; i < numpanels; i++) {
	 for (j = 0; j < charge->numcharges; j++) {
		Vector3D_sub(rvector, panels[i]->centroid, charge->points[j]);
		r = Vector3D_length(rvector);
		rhs[i] += charge->charges[j] * 1.0 / (4.0 * M_PI * idiel * r);
	 }
  }
  Vector3D_free(rvector);
}

void Charge_makerhs_ecf_qual(Vector rhs, Charge charge, Panel* panels, unsigned int numpanels, real idiel, real odiel) {
   unsigned int i, j;

   for (i = 0; i < numpanels; i++) {
      rhs[i] = (real)0;
      for (j = 0; j < charge->numcharges; j++) {
         real slp, dlp;

         dlp = Integration(charge->points[j], panels[i], POISSON_KERNEL, NULL, DOUBLE_LAYER_INT);
         rhs[i] += -dlp * charge->charges[j] / (4.0 * M_PI * idiel);
      }
  }

}

void Charge_makerhs_fromVector(Vector rhs, Charge charge, Vector sources, Panel* panels, unsigned int numpanels) {
   unsigned int i, j;
    
   for (i = 0; i < numpanels; i++) {
      rhs[i] = (real)0;
      for (j = 0; j < charge->numcharges; j++)
         rhs[i] += sources[j] / Vector3D_distance(charge->points[j], panels[i]->centroid);
   }
}
