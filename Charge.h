/* Typedef */

typedef struct _Charge {
  unsigned int numcharges;
  unsigned int globalindexstart;
  Vector3D* points;
  Vector charges;
} _Charge;

typedef _Charge* Charge;

/* Constructors and Destructors */

Charge Charge_allocate();
void Charge_free(Charge charge);

/* Operations */

void Charge_read(Charge charge, FILE* file);
void Charge_read_centroids(Charge charge, FILE* file, Vector3D* centroids, unsigned int numcentroids);

void Charge_makerhs(Vector rhs, Charge charge, Panel* panels, unsigned int numpanels);
void Charge_makerhs_ecf(Vector rhs, Charge charge, Panel* panels, unsigned int numpanels);
void Charge_makerhs_ecf_qual(Vector rhs, Charge charge, Panel* panels, unsigned int numpanels, real idiel, real odiel);
void Charge_makerhs_juffersimple(Vector rhs, Charge charge, Panel* panels, unsigned int numpanels, real idiel, real odiel);
void Charge_makerhs_fromVector(Vector rhs, Charge charge, Vector sources, Panel* panels, unsigned int numpanels);
