#include "FFTSVD.h"

/* Constructors and Destructors */

LJparameters LJparameters_allocate() {
   return (LJparameters)calloc(1, sizeof(_LJparameters));
}

void LJparameters_free(LJparameters ljparameters) {
   unsigned int i;
   
   for (i = 0; i < ljparameters->numatoms; i++)
      Vector3D_free(ljparameters->coordinates[i]);
   
   free(ljparameters->coordinates);
   Vector_free(ljparameters->R12terms);
   Vector_free(ljparameters->R6terms);

   free(ljparameters);
}

/* Operations */

void LJparameters_read(LJparameters ljparameters, FILE* file) {
   char line[128];
   unsigned int i;
        
   ljparameters->numatoms = 0;
    
   while (fgets(line, 128, file))
      ljparameters->numatoms++;

   ljparameters->coordinates = (Vector3D*)calloc(ljparameters->numatoms, sizeof(Vector3D));
   ljparameters->R12terms = Vector_allocate(ljparameters->numatoms);
   ljparameters->R6terms = Vector_allocate(ljparameters->numatoms);
    
   rewind(file);
    
   for (i = 0; i < ljparameters->numatoms; i++) {
      ljparameters->coordinates[i] = Vector3D_allocate();

      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf %lf %lf", &ljparameters->coordinates[i]->x, &ljparameters->coordinates[i]->y, &ljparameters->coordinates[i]->z, &ljparameters->R12terms[i], &ljparameters->R6terms[i]);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f %f %f", &ljparameters->coordinates[i]->x, &ljparameters->coordinates[i]->y, &ljparameters->coordinates[i]->z, &ljparameters->R12terms[i], &ljparameters->R6terms[i]);
#endif
#endif
   }
}

void LJparameters_makerhs(Vector rhs, LJparameters ljparameters, Panel* panels, unsigned int numpanels) {
   unsigned int i;
    
   for (i = 0; i < numpanels; i++) {
      rhs[i] = (real)1.0;
   }
}
