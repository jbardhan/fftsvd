/* Typedef */

typedef struct _QUI {
   unsigned int numpanels;
   unsigned int numconductors;
   unsigned int **conductorpanelindices;
   unsigned int *conductorpanelcounts;
   Vector3D* v1;
   Vector3D* v2;
   Vector3D* v3;
} _QUI;

typedef _QUI* QUI;

/* Constructors and Destructors */

QUI QUI_allocate();
void QUI_free(QUI qui);

/* Operations */

void QUI_read(QUI qui, FILE* file);

void QUI_getpanels(QUI qui, Panel* panels);
