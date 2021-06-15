#include "FFTSVD.h"

/* Constructors and Destructors */

QUI QUI_allocate() {
   return (QUI)calloc(1, sizeof(_QUI));
}

void QUI_free(QUI qui) {
   unsigned int i;
   
   for (i = 0; i < qui->numpanels; i++) {
      Vector3D_free(qui->v1[i]);
      Vector3D_free(qui->v2[i]);
      Vector3D_free(qui->v3[i]);
   }
   
   free(qui->conductorpanelindices);
   free(qui->v1);
   free(qui->v2);
   free(qui->v3);

   free(qui);
}

/* Operations */

void QUI_read(QUI qui, FILE* file) {
   char line[256];
   unsigned int i = 0, conductorindex;
   unsigned int *currentIndex;
   char c;
      
   qui->numpanels = 0;
   qui->numconductors = 0;
   while (fgets(line, 256, file)) {
      if ((line[0] == 'T') || (line[0] == 't')) {
         qui->numpanels++;
         sscanf(line, "%c %u", &c, &conductorindex);
         if (conductorindex > qui->numconductors)
            qui->numconductors = conductorindex;
      }
      else if ((line[0] == 'Q') || (line[0] == 'q')) {
         qui->numpanels += 2;
         sscanf(line, "%c %u", &c, &conductorindex);
         if (conductorindex > qui->numconductors)
            qui->numconductors = conductorindex;
      }
   }
   qui->conductorpanelindices = (unsigned int **)calloc(qui->numconductors, sizeof(unsigned int *));
   qui->conductorpanelcounts = (unsigned int*)calloc(qui->numconductors, sizeof(unsigned int));
   currentIndex = (unsigned int *) calloc(qui->numconductors, sizeof(unsigned int));
   rewind(file);

   while(fgets(line, 256, file)) {
      sscanf(line, "%c %u", &c, &conductorindex);  /* remember that the woven examples start with index 1!!!! */
      if (c == 'T' || c == 't')
         qui->conductorpanelcounts[conductorindex-1]++;
      else if (c == 'Q' || c == 'q')
         qui->conductorpanelcounts[conductorindex-1] += 2;
   }
   for (i=0; i < qui->numconductors; i++)
      qui->conductorpanelindices[i] = (unsigned int *)calloc(qui->conductorpanelcounts[i], sizeof(unsigned int));

   qui->v1 = (Vector3D*)calloc(qui->numpanels, sizeof(Vector3D));
   qui->v2 = (Vector3D*)calloc(qui->numpanels, sizeof(Vector3D));
   qui->v3 = (Vector3D*)calloc(qui->numpanels, sizeof(Vector3D));
   
   rewind(file);
   i = 0;
   while (fgets(line, 256, file)) {
      if ((line[0] == 'T') || (line[0] == 't')) {
         qui->v1[i] = Vector3D_allocate();
         qui->v2[i] = Vector3D_allocate();
         qui->v3[i] = Vector3D_allocate();

#ifdef REAL_IS_DOUBLE
         sscanf(line, "%c %u %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &c, &conductorindex,
                &qui->v1[i]->x, &qui->v1[i]->y, &qui->v1[i]->z,
                &qui->v2[i]->x, &qui->v2[i]->y, &qui->v2[i]->z,
                &qui->v3[i]->x, &qui->v3[i]->y, &qui->v3[i]->z);
#else
#ifdef REAL_IS_FLOAT
         sscanf(line, "%c %u %f %f %f %f %f %f %f %f %f",
                &c, &conductorindex,
                &qui->v1[i]->x, &qui->v1[i]->y, &qui->v1[i]->z,
                &qui->v2[i]->x, &qui->v2[i]->y, &qui->v2[i]->z,
                &qui->v3[i]->x, &qui->v3[i]->y, &qui->v3[i]->z);
#endif
#endif
         qui->conductorpanelindices[conductorindex-1][currentIndex[conductorindex-1]] = i;
         currentIndex[conductorindex-1]++;
         i++;
      }
      else if ((line[0] == 'Q') || (line[0] == 'q')) {
         Vector3D v4 = Vector3D_allocate();
         qui->v1[i] = Vector3D_allocate();
         qui->v2[i] = Vector3D_allocate();
         qui->v3[i] = Vector3D_allocate();
         qui->v1[i+1] = Vector3D_allocate();
         qui->v2[i+1] = Vector3D_allocate();
         qui->v3[i+1] = Vector3D_allocate();

#ifdef REAL_IS_DOUBLE
         sscanf(line, "%c %u %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &c, &conductorindex,
                &qui->v1[i]->x, &qui->v1[i]->y, &qui->v1[i]->z,
                &qui->v2[i]->x, &qui->v2[i]->y, &qui->v2[i]->z,
                &qui->v3[i]->x, &qui->v3[i]->y, &qui->v3[i]->z,
                &v4->x, &v4->y, &v4->z);
#else
#ifdef REAL_IS_FLOAT
         sscanf(line, "%c %u %f %f %f %f %f %f %f %f %f %f %f %f",
                &c, &conductorindex,
                &qui->v1[i]->x, &qui->v1[i]->y, &qui->v1[i]->z,
                &qui->v2[i]->x, &qui->v2[i]->y, &qui->v2[i]->z,
                &qui->v3[i]->x, &qui->v3[i]->y, &qui->v3[i]->z,
                &v4->x, &v4->y, &v4->z);
#endif
#endif
         Vector3D_copy(qui->v1[i+1], qui->v3[i]);
         Vector3D_copy(qui->v2[i+1], v4);
         Vector3D_copy(qui->v3[i+1], qui->v1[i]);

         Vector3D_free(v4);
         qui->conductorpanelindices[conductorindex-1][currentIndex[conductorindex-1]] = i;
         qui->conductorpanelindices[conductorindex-1][currentIndex[conductorindex-1]+1] = i+1;
         currentIndex[conductorindex-1]+=2;

         i += 2;
      }
   }
}

void QUI_getpanels(QUI qui, Panel* panels) {
   unsigned int i;
   
   for (i = 0; i < qui->numpanels; i++) {
      FlatPanel flatpanel = FlatPanel_allocate(qui->v1[i], qui->v2[i], qui->v3[i]);
      panels[i] = Panel_allocate();
      Panel_FlatPanel(panels[i], flatpanel);
   }
}
