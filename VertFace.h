/* Typedef */

typedef struct _VertFace {
   unsigned int numvertices, numfaces;
   Vector3D* vertices;
   unsigned int* facesx;
   unsigned int* facesy;
   unsigned int* facesz;
   unsigned int* facesgenus;
} _VertFace;

typedef _VertFace* VertFace;

/* Constructors and Destructors */

VertFace VertFace_allocate();
void VertFace_free(VertFace vf);

/* Operations */

void VertFace_readvert(VertFace vf, FILE* file);
void VertFace_readface(VertFace vf, FILE* file);
void VertFace_readface_flip(VertFace vf, FILE* file);

void VertFace_fix(VertFace vf, unsigned int accessible);

void VertFace_getpanels(VertFace vf, Panel* panels);
