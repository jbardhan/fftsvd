/* Typedef */

typedef struct _LJparameters {
   unsigned int numatoms;
   Vector3D* coordinates;
   Vector R12terms;
   Vector R6terms;
} _LJparameters;

typedef _LJparameters* LJparameters;

/* Constructors and Destructors */

LJparameters LJparameters_allocate();
void LJparameters_free(LJparameters ljparameters);

/* Operations */

void LJparameters_read(LJparameters ljparameters, FILE* file);

void LJparameters_makerhs(Vector rhs, LJparameters ljparameters, Panel* panels, unsigned int numpanels);
