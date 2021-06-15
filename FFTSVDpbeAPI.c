#include "FFTSVDpbeAPI.h"

/* Important Filenames */

const char* delphiparameterfilename = "fort.10";
const char* delphisizfilename = "fort.11";
const char* delphicrgfilename = "fort.12";
const char* delphipdbinputfilename = "fort.13";
const char* delphiphioutputfilename = "fort.14";
const char* delphipdbcollectfilename = "fort.15";
const char* delphifrcfilename = "fort.16";
const char* delphiphiinputfilename = "fort.18";
const char* delphipdboutputfilename = "fort.19";
const char* delphisrffilename = "fort.20";

/* Parameters */

real dielectricvertexdensity = 0.0f;
real saltvertexdensity = 0.0f;
real innerdielectric = 0.0f;
real outerdielectric = 0.0f;
real ionicstrength = 0.0f;
real ionexclusionradius = 0.0f;
real proberadius = 0.0f;
real errortolerance = 0.0f;
char energycalculation = '\0';
unsigned int outputmodpdb = 0;
unsigned int outputfrc = 0;
unsigned int addcoulombic = 0;
unsigned int accelerateM3 = 0;
unsigned int usequalocation = 0;
unsigned int usecurved = 0;
real inftydielectric = 0.0f;
real lambda          = 0.0f;

/* Utility */

void error(char* message, ...) {
   va_list ap;

   va_start(ap, message);

   printf("ERROR: ");
   vprintf(message, ap);
   printf("\n");

   va_end(ap);

   exit(-1);
}

void warning(char* message, ...) {
   va_list ap;

   va_start(ap, message);

   printf("WARNING: ");
   vprintf(message, ap);
   printf("\n");

   va_end(ap);
}

void removeWhitespace(char* s) {
   int length = strlen(s);
   char* temp = (char*)calloc(length+1, sizeof(char));
   int scount = 0, tempcount = 0;

   memset(temp, 0, length+1);

   for (scount = 0; scount < length; scount++)
      if (!isspace(s[scount])) {
         temp[tempcount] = s[scount];
         tempcount++;
      }

   strcpy(s, temp);

   free(temp);
}

int yesno(char* s) {
   if (strlen(s) == 0)
      return 0;
   else if (!strcasecmp(s, "t") || !strcasecmp(s, "y") || !strcasecmp(s, "true") || !strcasecmp(s, "yes") || !strcasecmp(s, "1"))
      return 1;
   else if (!strcasecmp(s, "f") || !strcasecmp(s, "n") || !strcasecmp(s, "false") || !strcasecmp(s, "no") || !strcasecmp(s, "0"))
      return 0;
   else {
      error("Invalid parameter boolean expression");
      return 0;
   }
}

void readParams(const char* filename) {
   FILE* parameterfile = NULL;
   char line[81], command[25], subcommand[10], subcharacter;
   int commandcount;
   char* nextposition;

   parameterfile = fopen(filename, "r");

   if (!parameterfile)
      error("Could not open parameter file %s: %s", filename, strerror(errno));

   while (fgets(line, 81, parameterfile)) {
      removeWhitespace(line);

      memset(command, 0, 25);
      commandcount = 0;

      for (commandcount = 0; isalpha(line[commandcount]); commandcount++)
         command[commandcount] = toupper(line[commandcount]);

      /* command now contains the directive to be processed */

      /* First, process simple equalities */

      nextposition = line + commandcount + 1;  /* After the = */

      if (!strcmp(command, "DVERT"))
         dielectricvertexdensity = atof(nextposition);
      else if (!strcmp(command, "SVERT"))
         saltvertexdensity = atof(nextposition);
      else if (!strcmp(command, "INDI"))
         innerdielectric = atof(nextposition);
      else if (!strcmp(command, "EXDI"))
         outerdielectric = atof(nextposition);
		else if (!strcmp(command, "INFTYDI"))
		  inftydielectric = atof(nextposition);
		else if (!strcmp(command, "LAMBDA"))
		  lambda = atof(nextposition);
      else if (!strcmp(command, "SALT"))
         ionicstrength = atof(nextposition);
      else if (!strcmp(command, "IONRAD"))
         ionexclusionradius = atof(nextposition);
      else if (!strcmp(command, "PRBRAD"))
         proberadius = atof(nextposition);
      else if (!strcmp(command, "DE"))
         errortolerance = atof(nextposition);
      else if (!strcmp(command, "ADDCOUL")) {
         if (yesno(nextposition))
            addcoulombic = 1;
      }
      else if (!strcmp(command, "USEQUAL")) {
         if (yesno(nextposition))
            usequalocation = 1;
      }
      else if (!strcmp(command, "CURVED")) {
         if (yesno(nextposition))
            usecurved = 1;
      }

      /* Next, process parsed directives */

      else if (!strcmp(command, "ENERGY")) {
         sscanf(nextposition-1, "(%c)", &subcharacter);
         energycalculation = toupper(subcharacter);
      }
      else if (!strcmp(command, "OUT")) {
         sscanf(nextposition-1, "(%[^)])", subcommand);
         if (!strcasecmp(subcommand, "MODPDB"))
            outputmodpdb = TRUE;
         else if (!strcasecmp(subcommand, "FRC"))
            outputfrc = TRUE;
         else if (!strcasecmp(subcommand, "PHI"));
         else
            warning("Unknown OUT parameter subcommand %s", subcommand);
      }
      else
         warning("Unknown parameter command %s", command);
   }

   fclose(parameterfile);
}

void checkParams() {
   if (dielectricvertexdensity < 0.0f)
      error("Dielectric vertex density must be > 0.0");
   if (saltvertexdensity < 0.0f)
      error("Salt vertex density must be > 0.0");
   if (innerdielectric < 1.0f)
      error("Inner dielectric must be >= 1.0");
   if (outerdielectric < 1.0f)
      error("Outer dielectric must be >= 1.0");
   if (innerdielectric > outerdielectric)
      warning("Inner dielectric is greater than outer dielectric");
   if (ionicstrength < 0.0f)
      error("Ionic strength must be >= 0.0");
   if (ionexclusionradius < 0.0f)
      error("Ionic exclusion radius must be >= 0.0");
   if (proberadius < 0.0f)
      error("Probe radius must be >= 0.0");
   if (errortolerance <= 0.0f)
      error("Error tolerance must be >= 0.0");
   if (energycalculation != 'G')
      error("Only energy calculation 'G' is currently supported");
   if (usequalocation && (ionicstrength > 0.0f))
      error("Qualocation does not work with ionic strength > 0.0");
}

void readPDB(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries) {
   FILE* pdbfile = NULL;
   char line[81], number[9];
   int atomcount = 0;

   pdbfile = fopen(filename, "r");

   if (!pdbfile)
      error("Could not open PDB file %s: %s", filename, strerror(errno));

   /* First, read through the file counting the number of ATOM or HETATM lines */

   while (fgets(line, 81, pdbfile)) {
      if (!strncasecmp(line, "ATOM", 4) || !strncasecmp(line, "HETATM", 6))
         atomcount++;
   }

   /* Allocate the memory */

   *PDBentries = (PDBentry*)calloc(atomcount, sizeof(PDBentry));

   if (!(*PDBentries))
      error("Error in memory allocation: %s", strerror(errno));

   /* Rewind the PDB file */

   rewind(pdbfile);

   /* Now read through again and store the entries */

   *numPDBentries = atomcount;

   atomcount = 0;

   while (fgets(line, 81, pdbfile)) {
      if (!strncasecmp(line, "ATOM", 4) || !strncasecmp(line, "HETATM", 6)) {
         strncpy((*PDBentries)[atomcount].record, line, 6);
         strncpy(number, line + 6, 5);
         number[5] = '\0';
         (*PDBentries)[atomcount].atomnumber = atoi(number);
         strncpy((*PDBentries)[atomcount].atomname, line + 12, 4);
         (*PDBentries)[atomcount].alternatelocation = line[16];
         strncpy((*PDBentries)[atomcount].residuename, line + 17, 3);
         (*PDBentries)[atomcount].chain = line[21];
         strncpy(number, line + 22, 4);
         number[4] = '\0';
         (*PDBentries)[atomcount].residuenumber = atoi(number);
         (*PDBentries)[atomcount].residueinsertion = line[26];
         strncpy(number, line + 30, 8);
         number[8] = '\0';
         (*PDBentries)[atomcount].x = atof(number);
         strncpy(number, line + 38, 8);
         number[8] = '\0';
         (*PDBentries)[atomcount].y = atof(number);
         strncpy(number, line + 46, 8);
         number[8] = '\0';
         (*PDBentries)[atomcount].z = atof(number);
         strncpy(number, line + 54, 6);
         number[6] = '\0';
         (*PDBentries)[atomcount].occupancy = atof(number);
         strncpy(number, line + 60, 6);
         number[6] = '\0';
         (*PDBentries)[atomcount].temperature = atof(number);
         strncpy(number, line + 67, 3);
         number[3] = '\0';
         (*PDBentries)[atomcount].footnotenumber = atof(number);

         removeWhitespace((*PDBentries)[atomcount].record);
         removeWhitespace((*PDBentries)[atomcount].atomname);
         removeWhitespace((*PDBentries)[atomcount].residuename);

         atomcount++;
      }
   }

   fclose(pdbfile);
}

void readCRD(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries) {
   FILE* crdfile = NULL;
   char line[81], number[11];
   int atomcount = 0;

   crdfile = fopen(filename, "r");

   if (!crdfile)
      error("Could not open CRD file %s: %s", filename, strerror(errno));

   /* First, read through the header lines */

   while (fgets(line, 81, crdfile))
      if (line[0] != '*')
         break;

   atomcount = atoi(line);

   /* Allocate the memory */

   *PDBentries = (PDBentry*)calloc(atomcount, sizeof(PDBentry));

   if (!(*PDBentries))
      error("Error in memory allocation: %s", strerror(errno));

   /* Now read through again and store the entries */

   *numPDBentries = atomcount;

   atomcount = 0;

   while (fgets(line, 81, crdfile)) {
      strncpy(number, line, 5);
      number[5] = '\0';
      (*PDBentries)[atomcount].atomnumber = atoi(number);
      strncpy(number, line + 5, 5);
      number[5] = '\0';
      (*PDBentries)[atomcount].residuenumber = atoi(number);
      strncpy((*PDBentries)[atomcount].residuename, line + 11, 3);
      strncpy((*PDBentries)[atomcount].atomname, line + 16, 4);
      strncpy(number, line + 20, 10);
      number[10] = '\0';
      (*PDBentries)[atomcount].x = atof(number);
      strncpy(number, line + 30, 10);
      number[10] = '\0';
      (*PDBentries)[atomcount].y = atof(number);
      strncpy(number, line + 40, 10);
      number[10] = '\0';
      (*PDBentries)[atomcount].z = atof(number);
      (*PDBentries)[atomcount].chain = line[51];
      strncpy(number, line + 61, 10);
      number[10] = '\0';
      (*PDBentries)[atomcount].temperature = atof(number);

      removeWhitespace((*PDBentries)[atomcount].atomname);
      removeWhitespace((*PDBentries)[atomcount].residuename);

      atomcount++;
   }

   fclose(crdfile);
}

void readXYZR(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries) {
   FILE* xyzrfile = NULL;
   char line[81], number[11];
   int atomcount = 0;

   xyzrfile = fopen(filename, "r");

   if (!xyzrfile)
      error("Could not open XYZR file %s: %s", filename, strerror(errno));

   /* First, read through the file counting the number of lines */

   while (fgets(line, 81, xyzrfile))
      atomcount++;

   /* Allocate the memory */

   *PDBentries = (PDBentry*)calloc(atomcount, sizeof(PDBentry));

   if (!(*PDBentries))
      error("Error in memory allocation: %s", strerror(errno));

   /* Rewind the PDB file */

   rewind(xyzrfile);

   /* Now read through again and store the entries */

   *numPDBentries = atomcount;

   atomcount = 0;

   while (fgets(line, 81, xyzrfile)) {
      sscanf(line, "%lf %lf %lf %lf", &(*PDBentries)[atomcount].x, &(*PDBentries)[atomcount].y, &(*PDBentries)[atomcount].z, &(*PDBentries)[atomcount].radius);
      atomcount++;
   }

   fclose(xyzrfile);
}

void readSIZ(const char* filename, unsigned int* numSIZentries, SIZentry** SIZentries) {
   FILE* sizfile = NULL;
   char line[81], number[9];
   int radiicount = 0;

   sizfile = fopen(filename, "r");

   if (!sizfile)
      error("Could not open SIZ file %s: %s", filename, strerror(errno));

   /* First, read through the file counting the number of non comment lines */

   fgets(line, 81, sizfile);  /* first line is junk */

   while (fgets(line, 81, sizfile)) {
      if ((line[0] != '!') && (strlen(line) > 0))
         radiicount++;
   }

   /* Allocate the memory */

   *SIZentries = (SIZentry*)calloc(radiicount, sizeof(SIZentry));

   if (!(*SIZentries))
      error("Error in memory allocation: %s", strerror(errno));

   /* Rewind the SIZ file */

   rewind(sizfile);

   /* Now read through again and store the entries */

   *numSIZentries = radiicount;

   radiicount = 0;

   fgets(line, 81, sizfile);  /* first line is junk */

   while (fgets(line, 81, sizfile)) {
      if ((line[0] != '!') && (strlen(line) > 0)) {
         strncpy((*SIZentries)[radiicount].atomlabel, line, 6);
         strncpy((*SIZentries)[radiicount].residuelabel, line + 6, 3);
         strncpy(number, line + 9, 8);
         number[8] = '\0';
         (*SIZentries)[radiicount].radius = atof(number);

         removeWhitespace((*SIZentries)[radiicount].atomlabel);
         removeWhitespace((*SIZentries)[radiicount].residuelabel);

         radiicount++;
      }
   }

   fclose(sizfile);
}

void readCRG(const char* filename, unsigned int* numCRGentries, CRGentry** CRGentries) {
#ifndef SCATTER
   FILE* crgfile = NULL;
   char line[81], number[9];
   int chargecount = 0;

   crgfile = fopen(filename, "r");

   if (!crgfile)
      error("Could not open CRG file %s: %s", filename, strerror(errno));

   /* First, read through the file counting the number of non comment lines */

   fgets(line, 81, crgfile);  /* first line is junk */

   while (fgets(line, 81, crgfile)) {
      if ((line[0] != '!') && (strlen(line) > 0))
         chargecount++;
   }

   /* Allocate the memory */
	if (!(*CRGentries))
	  free(*CRGentries);
	
   *CRGentries = (CRGentry*)calloc(chargecount, sizeof(CRGentry));

   if (!(*CRGentries))
      error("Error in memory allocation: %s", strerror(errno));

   /* Rewind the CRG file */

   rewind(crgfile);

   /* Now read through again and store the entries */

   *numCRGentries = chargecount;

   chargecount = 0;

   fgets(line, 81, crgfile);  /* first line is junk */

   while (fgets(line, 81, crgfile)) {
      if ((line[0] != '!') && (strlen(line) > 0)) {
         strncpy((*CRGentries)[chargecount].atomlabel, line, 6);
         strncpy((*CRGentries)[chargecount].residuelabel, line + 6, 3);
         strncpy(number, line + 9, 4);
         number[4] = '\0';
         (*CRGentries)[chargecount].residuenumber = atoi(number);
         (*CRGentries)[chargecount].chain = line[13];
         strncpy(number, line + 14, 8);
         number[8] = '\0';
         (*CRGentries)[chargecount].charge = atof(number);

         removeWhitespace((*CRGentries)[chargecount].atomlabel);
         removeWhitespace((*CRGentries)[chargecount].residuelabel);

         chargecount++;
      }
   }

   fclose(crgfile);
#else
	printf("Skipping readCRG because scattering calculation enabled!\n");
	*numCRGentries = 0;
	*CRGentries = NULL;
#endif
}

/* Output */

void writeMODPDB(const char* filename, PDBentry* PDBentries, unsigned int numPDBentries) {
   FILE* pdbfile = NULL;
   int a;

   pdbfile = fopen(filename, "w");

   if (!pdbfile)
      error("Could not open PDB file %s: %s", filename, strerror(errno));

   fprintf(pdbfile, "HEADER output from FFTSVDpbe\n");
   fprintf(pdbfile, "HEADER atom radii in columns 55-60\n");
   fprintf(pdbfile, "HEADER atom charges in columns 61-67\n");

   for (a = 0; a < numPDBentries; a++)
      fprintf(pdbfile, "%-6s%5d %4s%1c%3s %1c%4d%1c   %8.3f%8.3f%8.3f%6.2f%7.3f\n",
              PDBentries[a].record, PDBentries[a].atomnumber,
              PDBentries[a].atomname, PDBentries[a].alternatelocation,
              PDBentries[a].residuename, PDBentries[a].chain,
              PDBentries[a].residuenumber, PDBentries[a].residueinsertion,
              PDBentries[a].x, PDBentries[a].y,
              PDBentries[a].z, PDBentries[a].radius,
              PDBentries[a].charge);

   fclose(pdbfile);
}

void writeFRC(const char* filename, PDBentry* PDBentries, unsigned int numPDBentries) {
   FILE* frcfile = NULL;
   int a;

   frcfile = fopen(filename, "w");

   if (!frcfile)
      error("Could not open FRC file %s: %s", filename, strerror(errno));

   fprintf(frcfile, "FFTSVDpbe SITE POTENTIAL FILE\n");
   fprintf(frcfile, "grid size,percent fill: %3d %8f\n", 0, 0.0);
   fprintf(frcfile, "inner,outer dielectric: %8f %8f\n", innerdielectric, outerdielectric);
   fprintf(frcfile, "ionic strength (M): %8f\n", ionicstrength);
   fprintf(frcfile, "ion excl., probe radius: %8f %8f\n", ionexclusionradius, proberadius);
   fprintf(frcfile, "linear, nolinear iterations: %8d %8d\n", 0, 0);
   fprintf(frcfile, "boundary condition: %2d\n", 0);
   fprintf(frcfile, "Data Output: COORDINATES CHARGE POTENTIALS FIELDS\n");
   fprintf(frcfile, "title: FFTSVDpbe: BEM PBE Solver\n");
   fprintf(frcfile, "\n");
   fprintf(frcfile, "\n");
   fprintf(frcfile, "    ATOM COORDINATES (X,Y,Z)     CHARGE   GRID PT.    GRID FIELDS: (Ex, Ey, Ez)\n");
   for (a = 0; a < numPDBentries; a++)
      fprintf(frcfile, "%9.3f%9.3f%9.3f%8.4f%15.6f  %9.4f%9.4f%9.4f\n",
              PDBentries[a].x, PDBentries[a].y, PDBentries[a].z,
              PDBentries[a].charge, PDBentries[a].potential,
              0.0f, 0.0f, 0.0f);

   fprintf(frcfile, "total energy =   %6f  kt\n", 0.0);
   
   fclose(frcfile);
}

/* Molecule Handling */

void assignRadiiCharges(PDBentry* PDBentries, unsigned int numPDBentries, SIZentry* SIZentries, unsigned int numSIZentries, CRGentry* CRGentries, unsigned int numCRGentries) {
#ifndef SCATTER
   unsigned int a, r, c, matchlevel;

   for (a = 0; a < numPDBentries; a++) {
      matchlevel = 0;

      for (r = 0; r < numSIZentries; r++) {
         if ((!strcmp(PDBentries[a].atomname, SIZentries[r].atomlabel)) && (!strcmp(PDBentries[a].residuename, SIZentries[r].residuelabel))) {
            PDBentries[a].radius = SIZentries[r].radius;
            matchlevel = 3;
            break;  /* A complete match, so we can stop here */
         }
         else if ((strlen(SIZentries[r].residuelabel) == 0) && (!strncmp(PDBentries[a].atomname, SIZentries[r].atomlabel, strlen(SIZentries[r].atomlabel)))) {
            if (matchlevel < 2) {
               PDBentries[a].radius = SIZentries[r].radius;
               matchlevel = 2;
            }
         }
         else if ((strlen(SIZentries[r].residuelabel) == 0) && (PDBentries[a].atomname[0] == SIZentries[r].atomlabel[0]) && (strlen(SIZentries[r].atomlabel) == 1)) {
            if (matchlevel < 1) {
               PDBentries[a].radius = SIZentries[r].radius;
               matchlevel = 1;
            }
         }
      }

      if (matchlevel == 0) {
         warning("Atom %d %s %s %c has no radius entry, defaulting to 0.0",
                 PDBentries[a].atomnumber, PDBentries[a].atomname,
                 PDBentries[a].residuename, PDBentries[a].chain);

         PDBentries[a].radius = 0.0f;
      }
      else if (matchlevel == 1)
         warning("Atom %d %s %s %c using default element radius of %f",
                 PDBentries[a].atomnumber, PDBentries[a].atomname,
                 PDBentries[a].residuename, PDBentries[a].chain,
                 PDBentries[a].radius);

      matchlevel = 0;

      for (c = 0; c < numCRGentries; c++) {
         if ((CRGentries[c].chain == PDBentries[a].chain) && (CRGentries[c].residuenumber == PDBentries[a].residuenumber) && (!strcmp(PDBentries[a].residuename, CRGentries[c].residuelabel)) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            PDBentries[a].charge = CRGentries[c].charge;
            matchlevel = 8;
            break;  /* A complete match, so we can stop here */
         }
         else if ((CRGentries[c].residuenumber == PDBentries[a].residuenumber) && (!strcmp(PDBentries[a].residuename, CRGentries[c].residuelabel)) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 7) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 7;
            }
         }
         else if ((CRGentries[c].chain == PDBentries[a].chain) && (!strcmp(PDBentries[a].residuename, CRGentries[c].residuelabel)) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 6) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 6;
            }
         }
         else if ((!strcmp(PDBentries[a].residuename, CRGentries[c].residuelabel)) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 5) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 5;
            }
         }
         else if ((CRGentries[c].chain == PDBentries[a].chain) && (CRGentries[c].residuenumber == PDBentries[a].residuenumber) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 4) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 4;
            }
         }
         else if ((CRGentries[c].residuenumber == PDBentries[a].residuenumber) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 3) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 3;
            }
         }
         else if ((CRGentries[c].chain == PDBentries[a].chain) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 2) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 2;
            }
         }
         else if (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel)) {
            if (matchlevel < 1) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 1;
            }
         }
      }

      if ((matchlevel == 0) && (numCRGentries > 0)) {
         warning("Atom %d %s %s %c has no charge entry, defaulting to 0.0",
                 PDBentries[a].atomnumber, PDBentries[a].atomname,
                 PDBentries[a].residuename, PDBentries[a].chain);

         PDBentries[a].charge = 0.0f;
      }
   }
#else
	printf("Skipping assignRadiiCharges because scattering enabled.\n");
#endif 
}

/* Surfacing */

void getFlatPanels(const char* basefilename, Panel** panels, unsigned int* numpanels, unsigned int accessible, unsigned int cavity) {
   char vertfilename[1024], facefilename[1024];
   FILE* vertfile = NULL;
   FILE* facefile = NULL;
   VertFace vf;

   sprintf(vertfilename, "%s.vert", basefilename);
   sprintf(facefilename, "%s.face", basefilename);

   vertfile = fopen(vertfilename, "r");
   if (!vertfile)
      error("Could not open vert file %s: %s", vertfilename, strerror(errno));

   facefile = fopen(facefilename, "r");
   if (!facefile)
      error("Could not open face file %s: %s", facefilename, strerror(errno));

   vf = VertFace_allocate();

   VertFace_readvert(vf, vertfile);
   if (cavity && !usequalocation)
      VertFace_readface(vf, facefile);
   else
      VertFace_readface_flip(vf, facefile);

   VertFace_fix(vf, accessible);
	
   *numpanels = vf->numfaces;
   *panels = (Panel*)calloc(*numpanels, sizeof(Panel));

   VertFace_getpanels(vf, *panels);

   VertFace_free(vf);

   fclose(vertfile);
   fclose(facefile);
}

void getCurvedPanels(const char* basefilename, Panel** panels, unsigned int* numpanels, unsigned int accessible, unsigned int cavity) {
   char gstfilename[1024], torfilename[1024];
   FILE* gstfile = NULL;
   FILE* torfile = NULL;
   GST* gstpanels;
   TOR* torpanels;
   unsigned int numgstpanels = 0, numokgstpanels = 0, numtorpanels = 0, i, count = 0;

   sprintf(gstfilename, "%s.gst", basefilename);
   sprintf(torfilename, "%s.tor", basefilename);

   gstfile = fopen(gstfilename, "r");
   if (!gstfile)
      error("Could not open gst file %s: %s", gstfilename, strerror(errno));

   if (!accessible) {
      torfile = fopen(torfilename, "r");
      if (!torfile)
         error("Could not open tor file %s: %s", torfilename, strerror(errno));
   }

   GST_readfile(&numgstpanels, &gstpanels, gstfile, usequalocation);
   if (!accessible)
      TOR_readfile(&numtorpanels, &torpanels, torfile, usequalocation);

   for (i = 0; i < numgstpanels; i++) {
      Panel p = Panel_allocate();
      Panel_GST(p, gstpanels[i], 0);
      if (p->area > 1e-6)
         numokgstpanels++;
      Panel_free(p);
   }

   printf("GST panels with bad area: %u\n", numgstpanels-numokgstpanels);

   *numpanels = numokgstpanels + numtorpanels;

   *panels = (Panel*)calloc(*numpanels, sizeof(Panel));

   for (i = 0; i < numgstpanels; i++) {
      Panel p = Panel_allocate();
      Panel_GST(p, gstpanels[i], 0);
      if (p->area > 1e-6) {
         (*panels)[count] = p;
         count++;
      }
      else {
         printf("Remove GST %u with area %e\n", i+1, p->area);
         Panel_free(p);
      }
   }

   for (i = 0; i < numtorpanels; i++) {
      (*panels)[i+numokgstpanels] = Panel_allocate();
      Panel_TOR((*panels)[i+numokgstpanels], torpanels[i], 0);
   }

   fclose(gstfile);
   if (!accessible)
      fclose(torfile);
}

void readSRF(const char* filename, Panel*** saltpanels, unsigned int** numsaltpanels, unsigned int* numsalts,
Panel*** dielectricpanels, unsigned int** numdielectricpanels, unsigned int* numdielectrics, unsigned int** dielectricparent,
Panel*** dielectriccavitypanels, unsigned int** numdielectriccavitypanels, unsigned int* numdielectriccavities, unsigned int** dielectriccavityparent,
Panel*** saltcavitypanels, unsigned int** numsaltcavitypanels, unsigned int* numsaltcavities, unsigned int** saltcavityparent,
unsigned int* numtotalpanels) {
   FILE* srffile = fopen(filename, "r");
   char salttype, dielectrictype;
   char line[1024];
   // begin relative path load ( JPB 12/06/05 ): so that executable
   // can be run in a directory separate from geometry if the SRF file
   // contains ABSOLUTE PATHS, this change will make no difference.
   // if the SRF contains RELATIVE PATHS, then this change will allow
   // you to run your FFTSVDpbeAPI-based codes from any directory and
   // use a .SRF and assoc. geometry files from another directory
   char pathtosrf[1024];
   char unsafefilename[1024];
   char *ptrtoLastSlash;
   unsigned int numchar;
   // end relative path load 
   char* curfilename;
   char* curnum;
   char* buffer;
   unsigned int i, p;
   real area;

   if (!srffile)
      error("Could not open surface file %s: %s", filename, strerror(errno));

   // begin path modification stuff
   strcpy(unsafefilename, filename);
   ptrtoLastSlash = strrchr(unsafefilename, '/');
   numchar = 0;
   if (ptrtoLastSlash != NULL) 
      numchar = ptrtoLastSlash - unsafefilename + 1;
   unsafefilename[numchar] = '\0';

   // end path modification stuff

   fgets(line, 1024, srffile);
   sscanf(line, "%c", &salttype);
   fgets(line, 1024, srffile);
   sscanf(line, "%c", &dielectrictype);

   printf("Salt layers are %s panels\n", (salttype == 'f' ? "flat" : "curved"));
   printf("Dielectric layers are %s panels\n", (dielectrictype == 'f' ? "flat" : "curved"));

   *numtotalpanels = 0;

   /* Salts */

   fgets(line, 1024, srffile);
   buffer = strdup(line);

   *numsalts = 0;
   while ((curfilename = strsep(&buffer, " \n")))
      if (strlen(curfilename) > 0)
         (*numsalts)++;

   if (*numsalts > 0) {
      printf("Number of salt layers: %u\n", *numsalts);
      *saltpanels = (Panel**)calloc(*numsalts, sizeof(Panel*));
      *numsaltpanels = (unsigned int*)calloc(*numsalts, sizeof(unsigned int));
   }

   buffer = strdup(line);

   for (i = 0; i < *numsalts; i++) {
      curfilename = strsep(&buffer, " \n");
      // only modify filename if curfilename contains a RELATIVE PATH!
      if (curfilename[0] != '/') {
         pathtosrf[0] = '\0';
         strcat(pathtosrf, unsafefilename);
         strcat(pathtosrf, curfilename);
         curfilename = pathtosrf;
      }

      if (salttype == 'f')
         getFlatPanels(curfilename, &(*saltpanels)[i], &(*numsaltpanels)[i], 1, 0);
      else if (salttype == 'c')
         getCurvedPanels(curfilename, &(*saltpanels)[i], &(*numsaltpanels)[i], 1, 0);

      area = 0.0;
      for (p = 0; p < (*numsaltpanels)[i]; p++)
         area += (*saltpanels)[i][p]->area;

      printf("   Salt layer %u has %u panels (%f A^2)\n", i, (*numsaltpanels)[i], area);

      *numtotalpanels += (*numsaltpanels)[i];
   }

   /* Dielectrics */

   fgets(line, 1024, srffile);
   buffer = strdup(line);

   *numdielectrics = 0;
   while ((curfilename = strsep(&buffer, " \n")))
      if (strlen(curfilename) > 0)
         (*numdielectrics)++;

   if (*numdielectrics > 0) {
      printf("Number of dielectric layers: %u\n", *numdielectrics);
      *dielectricpanels = (Panel**)calloc(*numdielectrics, sizeof(Panel*));
      *numdielectricpanels = (unsigned int*)calloc(*numdielectrics, sizeof(unsigned int));
      *dielectricparent = (unsigned int*)calloc(*numdielectrics, sizeof(unsigned int));
   }

   buffer = strdup(line);

   for (i = 0; i < *numdielectrics; i++) {
      curfilename = strsep(&buffer, " \n");
      // only modify filename if curfilename contains a RELATIVE PATH!
      if (curfilename[0] != '/') {
         pathtosrf[0] = '\0';
         strcat(pathtosrf, unsafefilename);
         strcat(pathtosrf, curfilename);
         curfilename = pathtosrf;
      }

      if (dielectrictype == 'f')
         getFlatPanels(curfilename, &(*dielectricpanels)[i], &(*numdielectricpanels)[i], 0, 0);
      else if (dielectrictype == 'c')
         getCurvedPanels(curfilename, &(*dielectricpanels)[i], &(*numdielectricpanels)[i], 0, 0);

      area = 0.0;
      for (p = 0; p < (*numdielectricpanels)[i]; p++)
         area += (*dielectricpanels)[i][p]->area;

      printf("   Dielectric layer %u has %u panels (%f A^2)\n", i, (*numdielectricpanels)[i], area);

      *numtotalpanels += (*numdielectricpanels)[i];
   }

   /* Dielectric Cavities */

   fgets(line, 1024, srffile);
   buffer = strdup(line);

   *numdielectriccavities = 0;
   while ((curfilename = strsep(&buffer, " \n")))
      if (strlen(curfilename) > 0)
         (*numdielectriccavities)++;

   if (*numdielectriccavities > 0) {
      printf("Number of dielectric cavities: %u\n", *numdielectriccavities);
      *dielectriccavitypanels = (Panel**)calloc(*numdielectriccavities, sizeof(Panel*));
      *numdielectriccavitypanels = (unsigned int*)calloc(*numdielectriccavities, sizeof(unsigned int));
      *dielectriccavityparent = (unsigned int*)calloc(*numdielectriccavities, sizeof(unsigned int));
   }

   buffer = strdup(line);

   for (i = 0; i < *numdielectriccavities; i++) {
      curfilename = strsep(&buffer, " \n");
      // only modify filename if curfilename contains a RELATIVE PATH!
      if (curfilename[0] != '/') {
         pathtosrf[0] = '\0';
         strcat(pathtosrf, unsafefilename);
         strcat(pathtosrf, curfilename);
         curfilename = pathtosrf;
      }

      if (dielectrictype == 'f')
         getFlatPanels(curfilename, &(*dielectriccavitypanels)[i], &(*numdielectriccavitypanels)[i], 0, 1);
      else if (dielectrictype == 'c')
         getCurvedPanels(curfilename, &(*dielectriccavitypanels)[i], &(*numdielectriccavitypanels)[i], 0, 1);

      area = 0.0;
      for (p = 0; p < (*numdielectriccavitypanels)[i]; p++)
         area += (*dielectriccavitypanels)[i][p]->area;

      printf("   Dielectric cavity %u has %u panels (%f A^2)\n", i, (*numdielectriccavitypanels)[i], area);

      *numtotalpanels += (*numdielectriccavitypanels)[i];
   }

   /* Salt Cavities */

   fgets(line, 1024, srffile);
   buffer = strdup(line);

   *numsaltcavities = 0;
   while ((curfilename = strsep(&buffer, " \n")))
      if (strlen(curfilename) > 0)
         (*numsaltcavities)++;

   if (*numsaltcavities > 0) {
      printf("Number of salt cavities: %u\n", *numsaltcavities);
      *saltcavitypanels = (Panel**)calloc(*numsaltcavities, sizeof(Panel*));
      *numsaltcavitypanels = (unsigned int*)calloc(*numsaltcavities, sizeof(unsigned int));
      *saltcavityparent = (unsigned int*)calloc(*numsaltcavities, sizeof(unsigned int));
   }

   buffer = strdup(line);

   for (i = 0; i < *numsaltcavities; i++) {
      curfilename = strsep(&buffer, " \n");
      // only modify filename if curfilename contains a RELATIVE PATH!
      if (curfilename[0] != '/') {
         pathtosrf[0] = '\0';
         strcat(pathtosrf, unsafefilename);
         strcat(pathtosrf, curfilename);
         curfilename = pathtosrf;
      }

      if (salttype == 'f')
         getFlatPanels(curfilename, &(*saltcavitypanels)[i], &(*numsaltcavitypanels)[i], 1, 1);
      else if (salttype == 'c')
         getCurvedPanels(curfilename, &(*saltcavitypanels)[i], &(*numsaltcavitypanels)[i], 1, 1);

      area = 0.0;
      for (p = 0; p < (*numsaltcavitypanels)[i]; p++)
         area += (*saltcavitypanels)[i][p]->area;

      printf("   Salt cavity %u has %u panels (%f A^2)\n", i, (*numsaltcavitypanels)[i], area);

      *numtotalpanels += (*numsaltcavitypanels)[i];
   }

   /* Now do parents */

   fgets(line, 1024, srffile);
   buffer = strdup(line);
   for (i = 0; i < *numdielectrics; i++) {
      curnum = strsep(&buffer, " \n");
      (*dielectricparent)[i] = atoi(curnum);
      printf("Dielectric layer %u is inside salt layer %u\n", i, (*dielectricparent)[i]);
   }

   fgets(line, 1024, srffile);
   buffer = strdup(line);
   for (i = 0; i < *numdielectriccavities; i++) {
      curnum = strsep(&buffer, " \n");
      (*dielectriccavityparent)[i] = atoi(curnum);
      if ((*dielectriccavityparent)[i] == 999) {
         printf("Dielectric cavity %u marked as bad by meshmaker!  Removing from calculation\n", i);
         *numtotalpanels -= (*numdielectriccavitypanels)[i];
      }
      else
         printf("Dielectric cavity %u is inside dielectric layer %u\n", i, (*dielectriccavityparent)[i]);
   }

   fgets(line, 1024, srffile);
   buffer = strdup(line);
   for (i = 0; i < *numsaltcavities; i++) {
      curnum = strsep(&buffer, " \n");
      (*saltcavityparent)[i] = atoi(curnum);
		if ((*saltcavityparent)[i] == 999) { // new starts here
         printf("Salt cavity %u marked as bad by meshmaker!  Removing from calculation\n", i);
         *numtotalpanels -= (*numsaltcavitypanels)[i];
      }
      else // new ends here
        printf("Salt cavity %u is inside dielectric cavity %u\n", i, (*saltcavityparent)[i]);
   }

   fclose(srffile);
}

void generateFlatSRF(const char* pdbfilename, const char* sizfilename, const char* srffilename) {
   char command[1024];

   sprintf(command, "%s/meshmaker %s %s %s %f %f %f %f %u %u /tmp", PWD, pdbfilename, sizfilename, srffilename, proberadius, ionexclusionradius, dielectricvertexdensity, saltvertexdensity, 1, (ionicstrength > 0.0 && ionexclusionradius > 1e-6 ? 1 : 0));

   system(command);
}

void generateCurvedSRF(const char* pdbfilename, const char* sizfilename, const char* srffilename) {
   char command[1024];

   //sprintf(command, "%s/meshmaker %s %s %s %f %f %f %f %u %u /tmp", PWD, pdbfilename, sizfilename, srffilename, proberadius, ionexclusionradius, dielectricvertexdensity, saltvertexdensity, 4, (ionicstrength > 0.0 && ionexclusionradius > 1e-6 ? 4 : 0));
   sprintf(command, "%s/meshmaker %s %s %s %f %f %f %f %u %u .", PWD, pdbfilename, sizfilename, srffilename, proberadius, ionexclusionradius, dielectricvertexdensity, saltvertexdensity, 4, (ionicstrength > 0.0 && ionexclusionradius > 1e-6 ? 4 : 0));

   system(command);
}

/* Surface Operator */

void generateSurfaceOperator(SurfaceOperator* pbesurfaceoperator,
PDBentry* PDBentries, unsigned int numPDBentries,
Panel** saltpanels, unsigned int* numsaltpanels, unsigned int numsalts,
Panel** dielectricpanels, unsigned int* numdielectricpanels, unsigned int numdielectrics, unsigned int* dielectricparent,
Panel** dielectriccavitypanels, unsigned int* numdielectriccavitypanels, unsigned int numdielectriccavities, unsigned int* dielectriccavityparent,
Panel** saltcavitypanels, unsigned int* numsaltcavitypanels, unsigned int numsaltcavities, unsigned int* saltcavityparent,
unsigned int numtotalpanels) {
  real kappa = sqrt(ionicstrength) / DEBYE_CONSTANT;
  
  /* Stern case */

   if ((ionicstrength > 0.0) && (ionexclusionradius > 1e-6)) {
      printf("\nUsing Stern layer formulation\n");

      printf("\nConstructing operator OUTSIDE salts:\n");

      (*pbesurfaceoperator) = SurfaceOperator_allocate();
      (*pbesurfaceoperator)->kernel = HELMHOLTZ_KERNEL;
      (*pbesurfaceoperator)->epsilon = outerdielectric;
      (*pbesurfaceoperator)->kappa = kappa;
      (*pbesurfaceoperator)->mynumpanels = 0;
      (*pbesurfaceoperator)->numchildren = numsalts;
      (*pbesurfaceoperator)->children = (SurfaceOperator*)calloc(numsalts, sizeof(SurfaceOperator));
      (*pbesurfaceoperator)->parent = NULL;
      (*pbesurfaceoperator)->charges = NULL;

      unsigned int i, p, numtotalsaltpanels = 0;

      for (i = 0; i < numsalts; i++)
         numtotalsaltpanels += numsaltpanels[i];

      (*pbesurfaceoperator)->resultInternal = Vector_allocate(numtotalsaltpanels);
      (*pbesurfaceoperator)->resultExternal = Vector_allocate(numtotalsaltpanels);

      Panel* totalsaltpanels = (Panel*)calloc(numtotalsaltpanels, sizeof(Panel));
      Vector3D* totalsaltcentroids = (Vector3D*)calloc(numtotalsaltpanels, sizeof(Vector3D));

      unsigned int saltcount = 0;

      for (i = 0; i < numsalts; i++)
         for (p = 0; p < numsaltpanels[i]; p++) {
            totalsaltpanels[saltcount] = saltpanels[i][p];
            totalsaltcentroids[saltcount] = saltpanels[i][p]->centroid;
            saltcount++;
         }

      (*pbesurfaceoperator)->tree = Tree_allocate(totalsaltpanels, numtotalsaltpanels, totalsaltcentroids, numtotalsaltpanels, MAX_PANELS_PER_FINEST_CUBE, HELMHOLTZ_KERNEL, &kappa, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
      Tree_lists((*pbesurfaceoperator)->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
      Tree_fill((*pbesurfaceoperator)->tree);      
#endif
      Tree_memory((*pbesurfaceoperator)->tree);      

      for (i = 0; i < numsalts; i++) {
         printf("\nConstructing operator INSIDE salt %u, OUTSIDE dielectrics:\n", i);

         (*pbesurfaceoperator)->children[i] = SurfaceOperator_allocate();
         (*pbesurfaceoperator)->children[i]->kernel = POISSON_KERNEL;
         (*pbesurfaceoperator)->children[i]->epsilon = outerdielectric;
         (*pbesurfaceoperator)->children[i]->kappa = 0.0;
         (*pbesurfaceoperator)->children[i]->mynumpanels = numsaltpanels[i];
         (*pbesurfaceoperator)->children[i]->parent = (*pbesurfaceoperator);

         unsigned int j, numdielectricchildren = 0;
         unsigned int numtotalsaltanddielectricpanels = numsaltpanels[i];

         for (j = 0; j < numdielectrics; j++)
            if (dielectricparent[j] == i) {
               numdielectricchildren++;
#ifndef SCATTER
               numtotalsaltanddielectricpanels += numdielectricpanels[j];
#endif
            }

         if (numdielectricchildren > 0)
            (*pbesurfaceoperator)->children[i]->children = (SurfaceOperator*)calloc(numdielectricchildren, sizeof(SurfaceOperator));
         else
            (*pbesurfaceoperator)->children[i]->children = NULL;

         (*pbesurfaceoperator)->children[i]->numchildren = numdielectricchildren;

         (*pbesurfaceoperator)->children[i]->resultInternal = Vector_allocate(numtotalsaltanddielectricpanels);
         (*pbesurfaceoperator)->children[i]->resultExternal = Vector_allocate(numtotalsaltanddielectricpanels);

         Panel* totalsaltanddielectricpanels = (Panel*)calloc(numtotalsaltanddielectricpanels, sizeof(Panel));
         Vector3D* totalsaltanddielectriccentroids = (Vector3D*)calloc(numtotalsaltanddielectricpanels, sizeof(Vector3D));

         unsigned int panelcount = 0;

         for (p = 0; p < numsaltpanels[i]; p++) {
            totalsaltanddielectricpanels[panelcount] = saltpanels[i][p];
            totalsaltanddielectriccentroids[panelcount] = saltpanels[i][p]->centroid;
            panelcount++;
         }

#ifndef SCATTER
         for (j = 0; j < numdielectrics; j++)
            if (dielectricparent[j] == i)
               for (p = 0; p < numdielectricpanels[j]; p++) {
                  totalsaltanddielectricpanels[panelcount] = dielectricpanels[j][p];
                  totalsaltanddielectriccentroids[panelcount] = dielectricpanels[j][p]->centroid;
                  panelcount++;
               }
#endif
			
         (*pbesurfaceoperator)->children[i]->tree = Tree_allocate(totalsaltanddielectricpanels, numtotalsaltanddielectricpanels, totalsaltanddielectriccentroids, numtotalsaltanddielectricpanels, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
         Tree_lists((*pbesurfaceoperator)->children[i]->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
         Tree_fill((*pbesurfaceoperator)->children[i]->tree);      
#endif
         Tree_memory((*pbesurfaceoperator)->children[i]->tree);      

         unsigned int dielectricchildcount = 0;

         for (j = 0; j < numdielectrics; j++) {
            if (dielectricparent[j] != i)
               continue;

            printf("\nDielectric surface %u is inside salt surface %u\n", j, i);
            printf("Constructing operator INSIDE dielectric %u, OUTSIDE dielectric cavities\n", j);

            SurfaceOperator saltparent = (*pbesurfaceoperator)->children[i];

            saltparent->children[dielectricchildcount] = SurfaceOperator_allocate();
            saltparent->children[dielectricchildcount]->kernel = POISSON_KERNEL;
            saltparent->children[dielectricchildcount]->epsilon = innerdielectric;
            saltparent->children[dielectricchildcount]->kappa = 0.0;
            saltparent->children[dielectricchildcount]->mynumpanels = numdielectricpanels[j];
            saltparent->children[dielectricchildcount]->parent = saltparent;

            unsigned int k, numdielectriccavitychildren = 0;
            unsigned int numtotaldielectricandcavitypanels = numdielectricpanels[j];

            for (k = 0; k < numdielectriccavities; k++)
               if (dielectriccavityparent[k] == j) {
                  numdielectriccavitychildren++;
#ifndef SCATTER
                  numtotaldielectricandcavitypanels += numdielectriccavitypanels[k];
#endif
               }

            if (numdielectricchildren > 0)
               saltparent->children[dielectricchildcount]->children = (SurfaceOperator*)calloc(numdielectriccavitychildren, sizeof(SurfaceOperator));
            else
               saltparent->children[dielectricchildcount]->children = NULL;

            saltparent->children[dielectricchildcount]->numchildren = numdielectriccavitychildren;            

            saltparent->children[dielectricchildcount]->resultInternal = Vector_allocate(numdielectricpanels[j]);
            saltparent->children[dielectricchildcount]->resultExternal = Vector_allocate(numdielectricpanels[j]);

            Panel* totaldielectricandcavitypanels = (Panel*)calloc(numtotaldielectricandcavitypanels, sizeof(Panel));
            Vector3D* totaldielectricandcavitycentroids = (Vector3D*)calloc(numtotaldielectricandcavitypanels, sizeof(Vector3D));

            panelcount = 0;

            for (p = 0; p < numdielectricpanels[j]; p++) {
               totaldielectricandcavitypanels[panelcount] = dielectricpanels[j][p];
               totaldielectricandcavitycentroids[panelcount] = dielectricpanels[j][p]->centroid;
               panelcount++;
            }

#ifndef SCATTER
            for (k = 0; k < numdielectriccavities; k++)
               if (dielectriccavityparent[k] == j)
                  for (p = 0; p < numdielectriccavitypanels[k]; p++) {
                     totaldielectricandcavitypanels[panelcount] = dielectriccavitypanels[k][p];
                     totaldielectricandcavitycentroids[panelcount] = dielectriccavitypanels[k][p]->centroid;
                     panelcount++;
                  }
#endif
            saltparent->children[dielectricchildcount]->tree = Tree_allocate(totaldielectricandcavitypanels, numtotaldielectricandcavitypanels, totaldielectricandcavitycentroids, numtotaldielectricandcavitypanels, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
            Tree_lists(saltparent->children[dielectricchildcount]->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
            Tree_fill(saltparent->children[dielectricchildcount]->tree);
#endif
            Tree_memory(saltparent->children[dielectricchildcount]->tree);
            
            /* There might be charges in this layer, so handle them */

            unsigned int a, numcharges = numPDBentries;

            saltparent->children[dielectricchildcount]->charges = Charge_allocate();
            saltparent->children[dielectricchildcount]->charges->numcharges = numcharges;
            saltparent->children[dielectricchildcount]->charges->globalindexstart = 0;
            saltparent->children[dielectricchildcount]->charges->points = (Vector3D*)calloc(numcharges, sizeof(Vector3D));
            saltparent->children[dielectricchildcount]->charges->charges = Vector_allocate(numcharges);

            unsigned int chargecount = 0;

            for (a = 0; a < numPDBentries; a++) {
               Vector3D charge = Vector3D_allocate();
               charge->x = PDBentries[a].x;
               charge->y = PDBentries[a].y;
               charge->z = PDBentries[a].z;
               saltparent->children[dielectricchildcount]->charges->points[chargecount] = charge;
               if (toupper(PDBentries[a].chain) != 'X')
                  saltparent->children[dielectricchildcount]->charges->charges[chargecount] = PDBentries[a].charge;
               else
                  saltparent->children[dielectricchildcount]->charges->charges[chargecount] = 0.0;
               chargecount++;
            }

            // Create M3 tree if requested

            if (accelerateM3) {
               printf("\nConstructing M3 operator INSIDE dielectric %u, OUTSIDE dielectric cavities\n", j);
               saltparent->children[dielectricchildcount]->M3 = Tree_allocate(totaldielectricandcavitypanels, numtotaldielectricandcavitypanels, saltparent->children[dielectricchildcount]->charges->points, numcharges, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
               Tree_lists(saltparent->children[dielectricchildcount]->M3);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
               Tree_fill(saltparent->children[dielectricchildcount]->M3);
#endif
               Tree_memory(saltparent->children[dielectricchildcount]->M3);
            }

            unsigned int dielectriccavitychildcount = 0;

            for (k = 0; k < numdielectriccavities; k++) {
               if (dielectriccavityparent[k] != j)
                  continue;

               printf("\nDielectric cavity surface %u is inside dielectric surface %u\n", k, j);
               printf("Constructing operator INSIDE dielectric cavity %u, OUTSIDE salt cavities\n", k);

               SurfaceOperator dielectricparent = saltparent->children[dielectricchildcount];

               dielectricparent->children[dielectriccavitychildcount] = SurfaceOperator_allocate();
               dielectricparent->children[dielectriccavitychildcount]->kernel = POISSON_KERNEL;
               dielectricparent->children[dielectriccavitychildcount]->epsilon = outerdielectric;
               dielectricparent->children[dielectriccavitychildcount]->kappa = 0.0;
               dielectricparent->children[dielectriccavitychildcount]->mynumpanels = numdielectriccavitypanels[k];
               dielectricparent->children[dielectriccavitychildcount]->parent = dielectricparent;

               unsigned int l, numsaltcavitychildren = 0;
               unsigned int numtotaldielectricandsaltcavitypanels = numdielectriccavitypanels[k];

               for (l = 0; l < numsaltcavities; l++)
                  if (saltcavityparent[l] == k) {
                     numsaltcavitychildren++;
#ifndef SCATTER
                     numtotaldielectricandsaltcavitypanels += numsaltcavitypanels[l];
#endif
                  }

               if (numsaltcavitychildren > 0)
                  dielectricparent->children[dielectriccavitychildcount]->children = (SurfaceOperator*)calloc(numsaltcavitychildren, sizeof(SurfaceOperator));
               else
                  dielectricparent->children[dielectriccavitychildcount]->children = NULL;

               dielectricparent->children[dielectriccavitychildcount]->numchildren = numsaltcavitychildren;

               dielectricparent->children[dielectriccavitychildcount]->resultInternal = Vector_allocate(numdielectriccavitypanels[k]);
               dielectricparent->children[dielectriccavitychildcount]->resultExternal = Vector_allocate(numdielectriccavitypanels[k]);

               Panel* totaldielectricandsaltcavitypanels = (Panel*)calloc(numtotaldielectricandsaltcavitypanels, sizeof(Panel));
               Vector3D* totaldielectricandsaltcavitycentroids = (Vector3D*)calloc(numtotaldielectricandsaltcavitypanels, sizeof(Vector3D));

               panelcount = 0;

               for (p = 0; p < numdielectriccavitypanels[k]; p++) {
                  totaldielectricandsaltcavitypanels[panelcount] = dielectriccavitypanels[k][p];
                  totaldielectricandsaltcavitycentroids[panelcount] = dielectriccavitypanels[k][p]->centroid;
                  panelcount++;
               }

               for (l = 0; l < numsaltcavities; l++)
                  if (saltcavityparent[l] == k)
                     for (p = 0; p < numsaltcavitypanels[l]; p++) {
                        totaldielectricandsaltcavitypanels[panelcount] = saltcavitypanels[l][p];
                        totaldielectricandsaltcavitycentroids[panelcount] = saltcavitypanels[l][p]->centroid;                        
                        panelcount++;
                     }

               dielectricparent->children[dielectriccavitychildcount]->tree = Tree_allocate(totaldielectricandsaltcavitypanels, numtotaldielectricandsaltcavitypanels, totaldielectricandsaltcavitycentroids, numtotaldielectricandsaltcavitypanels, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
               Tree_lists(dielectricparent->children[dielectriccavitychildcount]->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
               Tree_fill(dielectricparent->children[dielectriccavitychildcount]->tree);
#endif
               Tree_memory(dielectricparent->children[dielectriccavitychildcount]->tree);

               unsigned int saltcavitychildcount = 0;

               for (l = 0; l < numsaltcavities; l++) {
                  if (saltcavityparent[l] != k)
                     continue;

                  printf("\nSalt cavity surface %u is inside dielectric cavity surface %u\n", l, k);
                  printf("Constructing operator INSIDE salt cavity %u\n", l);

                  SurfaceOperator dielectriccavityparent = dielectricparent->children[dielectriccavitychildcount];

                  dielectriccavityparent->children[saltcavitychildcount] = SurfaceOperator_allocate();
                  dielectriccavityparent->children[saltcavitychildcount]->kernel = HELMHOLTZ_KERNEL;
                  dielectriccavityparent->children[saltcavitychildcount]->epsilon = outerdielectric;
                  dielectriccavityparent->children[saltcavitychildcount]->kappa = kappa;
                  dielectriccavityparent->children[saltcavitychildcount]->mynumpanels = numsaltcavitypanels[l];
                  dielectriccavityparent->children[saltcavitychildcount]->numchildren = 0;
                  dielectriccavityparent->children[saltcavitychildcount]->children = NULL;
                  dielectriccavityparent->children[saltcavitychildcount]->parent = dielectriccavityparent;

                  dielectriccavityparent->children[saltcavitychildcount]->resultInternal = Vector_allocate(numsaltcavitypanels[l]);
                  dielectriccavityparent->children[saltcavitychildcount]->resultExternal = Vector_allocate(numsaltcavitypanels[l]);

                  Vector3D* saltcavitycentroids = (Vector3D*)calloc(numsaltcavitypanels[l], sizeof(Vector3D));

                  for (p = 0; p < numsaltcavitypanels[l]; p++)
                     saltcavitycentroids[p] = saltcavitypanels[l][p]->centroid;

                  dielectriccavityparent->children[saltcavitychildcount]->tree = Tree_allocate(saltcavitypanels[l], numsaltcavitypanels[l], saltcavitycentroids, numsaltcavitypanels[l], MAX_PANELS_PER_FINEST_CUBE, HELMHOLTZ_KERNEL, &kappa, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
                  Tree_lists(dielectriccavityparent->children[saltcavitychildcount]->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
                  Tree_fill(dielectriccavityparent->children[saltcavitychildcount]->tree);
#endif
                  Tree_memory(dielectriccavityparent->children[saltcavitychildcount]->tree);

                  saltcavitychildcount++;
               }

               dielectriccavitychildcount++;
            }

            dielectricchildcount++;
         }
      }
   }

   /* Salt, but no Stern */

   else if (ionicstrength > 0.0f) {
      printf("\nUsing Salt without Stern layer formulation\n");

      printf("\nConstructing operator OUTSIDE dielectrics:\n");

      (*pbesurfaceoperator) = SurfaceOperator_allocate();
      (*pbesurfaceoperator)->kernel = HELMHOLTZ_KERNEL;
      (*pbesurfaceoperator)->epsilon = outerdielectric;
      (*pbesurfaceoperator)->kappa = kappa;
      (*pbesurfaceoperator)->mynumpanels = 0;
      (*pbesurfaceoperator)->numchildren = numdielectrics;
      (*pbesurfaceoperator)->children = (SurfaceOperator*)calloc(numdielectrics, sizeof(SurfaceOperator));
      (*pbesurfaceoperator)->parent = NULL;
      (*pbesurfaceoperator)->charges = NULL;

      unsigned int i, p, numtotaldielectricpanels = 0;

      for (i = 0; i < numdielectrics; i++)
         numtotaldielectricpanels += numdielectricpanels[i];

      (*pbesurfaceoperator)->resultInternal = Vector_allocate(numtotaldielectricpanels);
      (*pbesurfaceoperator)->resultExternal = Vector_allocate(numtotaldielectricpanels);

      Panel* totaldielectricpanels = (Panel*)calloc(numtotaldielectricpanels, sizeof(Panel));
      Vector3D* totaldielectriccentroids = (Vector3D*)calloc(numtotaldielectricpanels, sizeof(Vector3D));

      unsigned int dielectriccount = 0;

      for (i = 0; i < numdielectrics; i++)
         for (p = 0; p < numdielectricpanels[i]; p++) {
            totaldielectricpanels[dielectriccount] = dielectricpanels[i][p];
            totaldielectriccentroids[dielectriccount] = dielectricpanels[i][p]->centroid;
            dielectriccount++;
         }

      (*pbesurfaceoperator)->tree = Tree_allocate(totaldielectricpanels, numtotaldielectricpanels, totaldielectriccentroids, numtotaldielectricpanels, MAX_PANELS_PER_FINEST_CUBE, HELMHOLTZ_KERNEL, &kappa, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
      Tree_lists((*pbesurfaceoperator)->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
      Tree_fill((*pbesurfaceoperator)->tree);      
#endif
      Tree_memory((*pbesurfaceoperator)->tree);      

      for (i = 0; i < numdielectrics; i++) {
         printf("\nConstructing operator INSIDE dielectric %u, OUTSIDE cavities:\n", i);

         (*pbesurfaceoperator)->children[i] = SurfaceOperator_allocate();
         (*pbesurfaceoperator)->children[i]->kernel = POISSON_KERNEL;
         (*pbesurfaceoperator)->children[i]->epsilon = innerdielectric;
         (*pbesurfaceoperator)->children[i]->kappa = 0.0;
         (*pbesurfaceoperator)->children[i]->mynumpanels = numdielectricpanels[i];
         (*pbesurfaceoperator)->children[i]->parent = (*pbesurfaceoperator);

         unsigned int j, numdielectriccavitychildren = 0;
         unsigned int numtotaldielectricandcavitypanels = numdielectricpanels[i];

         for (j = 0; j < numdielectriccavities; j++)
            if (dielectriccavityparent[j] == i) {
               numdielectriccavitychildren++;
               numtotaldielectricandcavitypanels += numdielectriccavitypanels[j];
            }

         if (numdielectriccavitychildren > 0)
            (*pbesurfaceoperator)->children[i]->children = (SurfaceOperator*)calloc(numdielectriccavitychildren, sizeof(SurfaceOperator));
         else
            (*pbesurfaceoperator)->children[i]->children = NULL;

         (*pbesurfaceoperator)->children[i]->numchildren = numdielectriccavitychildren;

         (*pbesurfaceoperator)->children[i]->resultInternal = Vector_allocate(numtotaldielectricandcavitypanels);
         (*pbesurfaceoperator)->children[i]->resultExternal = Vector_allocate(numtotaldielectricandcavitypanels);

         Panel* totaldielectricandcavitypanels = (Panel*)calloc(numtotaldielectricandcavitypanels, sizeof(Panel));
         Vector3D* totaldielectricandcavitycentroids = (Vector3D*)calloc(numtotaldielectricandcavitypanels, sizeof(Vector3D));

         unsigned int panelcount = 0;

         for (p = 0; p < numdielectricpanels[i]; p++) {
            totaldielectricandcavitypanels[panelcount] = dielectricpanels[i][p];
            totaldielectricandcavitycentroids[panelcount] = dielectricpanels[i][p]->centroid;
            panelcount++;
         }

         for (j = 0; j < numdielectriccavities; j++)
            if (dielectriccavityparent[j] == i)
               for (p = 0; p < numdielectriccavitypanels[j]; p++) {
                  totaldielectricandcavitypanels[panelcount] = dielectriccavitypanels[j][p];
                  totaldielectricandcavitycentroids[panelcount] = dielectriccavitypanels[j][p]->centroid;
                  panelcount++;
               }

         (*pbesurfaceoperator)->children[i]->tree = Tree_allocate(totaldielectricandcavitypanels, numtotaldielectricandcavitypanels, totaldielectricandcavitycentroids, numtotaldielectricandcavitypanels, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
         Tree_lists((*pbesurfaceoperator)->children[i]->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
         Tree_fill((*pbesurfaceoperator)->children[i]->tree);      
#endif
         Tree_memory((*pbesurfaceoperator)->children[i]->tree);      

         /* There might be charges in this layer, so handle them */

         unsigned int a, numcharges = numPDBentries;

         (*pbesurfaceoperator)->children[i]->charges = Charge_allocate();
         (*pbesurfaceoperator)->children[i]->charges->numcharges = numcharges;
         (*pbesurfaceoperator)->children[i]->charges->globalindexstart = 0;
         (*pbesurfaceoperator)->children[i]->charges->points = (Vector3D*)calloc(numcharges, sizeof(Vector3D));
         (*pbesurfaceoperator)->children[i]->charges->charges = Vector_allocate(numcharges);

         unsigned int chargecount = 0;

         for (a = 0; a < numPDBentries; a++) {
            Vector3D charge = Vector3D_allocate();
            charge->x = PDBentries[a].x;
            charge->y = PDBentries[a].y;
            charge->z = PDBentries[a].z;
            (*pbesurfaceoperator)->children[i]->charges->points[chargecount] = charge;
            if (toupper(PDBentries[a].chain) != 'X')
               (*pbesurfaceoperator)->children[i]->charges->charges[chargecount] = PDBentries[a].charge;
            else
               (*pbesurfaceoperator)->children[i]->charges->charges[chargecount] = 0.0;
            chargecount++;
         }

         // Create M3 tree if requested

         if (accelerateM3) {
            printf("\nConstructing M3 operator INSIDE dielectric %u, OUTSIDE dielectric cavities\n", i);
            (*pbesurfaceoperator)->children[i]->M3 = Tree_allocate(totaldielectricandcavitypanels, numtotaldielectricandcavitypanels, (*pbesurfaceoperator)->children[i]->charges->points, numcharges, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
            Tree_lists((*pbesurfaceoperator)->children[i]->M3);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
            Tree_fill((*pbesurfaceoperator)->children[i]->M3);
#endif
            Tree_memory((*pbesurfaceoperator)->children[i]->M3);
         }

         unsigned int dielectriccavitychildcount = 0;

         for (j = 0; j < numdielectriccavities; j++) {
            if (dielectriccavityparent[j] != i)
               continue;

            printf("\nDielectric cavity surface %u is inside dielectric surface %u\n", j, i);
            printf("Constructing operator INSIDE dielectric cavity %u\n", j);

            SurfaceOperator dielectricparent = (*pbesurfaceoperator)->children[i];

            dielectricparent->children[dielectriccavitychildcount] = SurfaceOperator_allocate();
            dielectricparent->children[dielectriccavitychildcount]->kernel = HELMHOLTZ_KERNEL;
            dielectricparent->children[dielectriccavitychildcount]->epsilon = outerdielectric;
            dielectricparent->children[dielectriccavitychildcount]->kappa = kappa;
            dielectricparent->children[dielectriccavitychildcount]->mynumpanels = numdielectriccavitypanels[j];
            dielectricparent->children[dielectriccavitychildcount]->numchildren = 0;
            dielectricparent->children[dielectriccavitychildcount]->children = NULL;
            dielectricparent->children[dielectriccavitychildcount]->parent = dielectricparent;
            dielectricparent->children[dielectriccavitychildcount]->charges = NULL;

            dielectricparent->children[dielectriccavitychildcount]->resultInternal = Vector_allocate(numdielectriccavitypanels[j]);
            dielectricparent->children[dielectriccavitychildcount]->resultExternal = Vector_allocate(numdielectriccavitypanels[j]);

            Vector3D* dielectriccavitycentroids = (Vector3D*)calloc(numdielectriccavitypanels[j], sizeof(Vector3D));

            for (p = 0; p < numdielectriccavitypanels[j]; p++)
               dielectriccavitycentroids[p] = dielectriccavitypanels[j][p]->centroid;

            dielectricparent->children[dielectriccavitychildcount]->tree = Tree_allocate(dielectriccavitypanels[j], numdielectriccavitypanels[j], dielectriccavitycentroids, numdielectriccavitypanels[j], MAX_PANELS_PER_FINEST_CUBE, HELMHOLTZ_KERNEL, &kappa, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
            Tree_lists(dielectricparent->children[dielectriccavitychildcount]->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
            Tree_fill(dielectricparent->children[dielectriccavitychildcount]->tree);
#endif
            Tree_memory(dielectricparent->children[dielectriccavitychildcount]->tree);

            dielectriccavitychildcount++;
         }
      }
   }

   /* No salt, just dielectrics */
   else {
      printf("\nUsing no Salt formulation\n");

      printf("\nConstructing operator OUTSIDE dielectrics:\n");

      (*pbesurfaceoperator) = SurfaceOperator_allocate();
      (*pbesurfaceoperator)->kernel = POISSON_KERNEL;
      (*pbesurfaceoperator)->epsilon = outerdielectric;
      (*pbesurfaceoperator)->kappa = 0.0;
      (*pbesurfaceoperator)->mynumpanels = 0;
      (*pbesurfaceoperator)->numchildren = numdielectrics;
      (*pbesurfaceoperator)->children = (SurfaceOperator*)calloc(numdielectrics, sizeof(SurfaceOperator));
      (*pbesurfaceoperator)->parent = NULL;
      (*pbesurfaceoperator)->charges = NULL;

      unsigned int i, p, numtotaldielectricpanels = 0;

      for (i = 0; i < numdielectrics; i++)
         numtotaldielectricpanels += numdielectricpanels[i];

      (*pbesurfaceoperator)->resultInternal = Vector_allocate(numtotaldielectricpanels);
      (*pbesurfaceoperator)->resultExternal = Vector_allocate(numtotaldielectricpanels);

      Panel* totaldielectricpanels = (Panel*)calloc(numtotaldielectricpanels, sizeof(Panel));
      Vector3D* totaldielectriccentroids = (Vector3D*)calloc(numtotaldielectricpanels, sizeof(Vector3D));

      unsigned int dielectriccount = 0;

      for (i = 0; i < numdielectrics; i++)
         for (p = 0; p < numdielectricpanels[i]; p++) {
            totaldielectricpanels[dielectriccount] = dielectricpanels[i][p];
            totaldielectriccentroids[dielectriccount] = dielectricpanels[i][p]->centroid;
            dielectriccount++;
         }

      (*pbesurfaceoperator)->tree = Tree_allocate(totaldielectricpanels, numtotaldielectricpanels, totaldielectriccentroids, numtotaldielectricpanels, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
      Tree_lists((*pbesurfaceoperator)->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
      Tree_fill((*pbesurfaceoperator)->tree);      
#endif
      Tree_memory((*pbesurfaceoperator)->tree);      

      for (i = 0; i < numdielectrics; i++) {
         printf("\nConstructing operator INSIDE dielectric %u, OUTSIDE cavities:\n", i);

         (*pbesurfaceoperator)->children[i] = SurfaceOperator_allocate();
         (*pbesurfaceoperator)->children[i]->kernel = POISSON_KERNEL;
         (*pbesurfaceoperator)->children[i]->epsilon = innerdielectric;
         (*pbesurfaceoperator)->children[i]->kappa = 0.0;
         (*pbesurfaceoperator)->children[i]->mynumpanels = numdielectricpanels[i];
         (*pbesurfaceoperator)->children[i]->parent = (*pbesurfaceoperator);

         unsigned int j, numdielectriccavitychildren = 0;
         unsigned int numtotaldielectricandcavitypanels = numdielectricpanels[i];

         for (j = 0; j < numdielectriccavities; j++)
            if (dielectriccavityparent[j] == i) {
               numdielectriccavitychildren++;
               numtotaldielectricandcavitypanels += numdielectriccavitypanels[j];
            }

         if (numdielectriccavitychildren > 0)
            (*pbesurfaceoperator)->children[i]->children = (SurfaceOperator*)calloc(numdielectriccavitychildren, sizeof(SurfaceOperator));
         else
            (*pbesurfaceoperator)->children[i]->children = NULL;

         (*pbesurfaceoperator)->children[i]->numchildren = numdielectriccavitychildren;

         (*pbesurfaceoperator)->children[i]->resultInternal = Vector_allocate(numtotaldielectricandcavitypanels);
         (*pbesurfaceoperator)->children[i]->resultExternal = Vector_allocate(numtotaldielectricandcavitypanels);

         Panel* totaldielectricandcavitypanels = (Panel*)calloc(numtotaldielectricandcavitypanels, sizeof(Panel));
         Vector3D* totaldielectricandcavitycentroids = (Vector3D*)calloc(numtotaldielectricandcavitypanels, sizeof(Vector3D));

         unsigned int panelcount = 0;

         for (p = 0; p < numdielectricpanels[i]; p++) {
            totaldielectricandcavitypanels[panelcount] = dielectricpanels[i][p];
            totaldielectricandcavitycentroids[panelcount] = dielectricpanels[i][p]->centroid;
            panelcount++;
         }

         for (j = 0; j < numdielectriccavities; j++)
            if (dielectriccavityparent[j] == i)
               for (p = 0; p < numdielectriccavitypanels[j]; p++) {
                  totaldielectricandcavitypanels[panelcount] = dielectriccavitypanels[j][p];
                  totaldielectricandcavitycentroids[panelcount] = dielectriccavitypanels[j][p]->centroid;
                  panelcount++;
               }

         (*pbesurfaceoperator)->children[i]->tree = Tree_allocate(totaldielectricandcavitypanels, numtotaldielectricandcavitypanels, totaldielectricandcavitycentroids, numtotaldielectricandcavitypanels, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
         Tree_lists((*pbesurfaceoperator)->children[i]->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
         Tree_fill((*pbesurfaceoperator)->children[i]->tree);      
#endif
         Tree_memory((*pbesurfaceoperator)->children[i]->tree);      

         /* There might be charges in this layer, so handle them */

         unsigned int a, numcharges = numPDBentries;

         (*pbesurfaceoperator)->children[i]->charges = Charge_allocate();
         (*pbesurfaceoperator)->children[i]->charges->numcharges = numcharges;
         (*pbesurfaceoperator)->children[i]->charges->globalindexstart = 0;
         (*pbesurfaceoperator)->children[i]->charges->points = (Vector3D*)calloc(numcharges, sizeof(Vector3D));
         (*pbesurfaceoperator)->children[i]->charges->charges = Vector_allocate(numcharges);

         unsigned int chargecount = 0;

         for (a = 0; a < numPDBentries; a++) {
            Vector3D charge = Vector3D_allocate();
            charge->x = PDBentries[a].x;
            charge->y = PDBentries[a].y;
            charge->z = PDBentries[a].z;
            (*pbesurfaceoperator)->children[i]->charges->points[chargecount] = charge;
            if (toupper(PDBentries[a].chain) != 'X')
               (*pbesurfaceoperator)->children[i]->charges->charges[chargecount] = PDBentries[a].charge;
            else
               (*pbesurfaceoperator)->children[i]->charges->charges[chargecount] = 0.0;
            chargecount++;
         }

         // Create M3 tree if requested

         if (accelerateM3) {
            printf("\nConstructing M3 operator INSIDE dielectric %u, OUTSIDE dielectric cavities\n", i);
            (*pbesurfaceoperator)->children[i]->M3 = Tree_allocate(totaldielectricandcavitypanels, numtotaldielectricandcavitypanels, (*pbesurfaceoperator)->children[i]->charges->points, numcharges, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
            Tree_lists((*pbesurfaceoperator)->children[i]->M3);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
            Tree_fill((*pbesurfaceoperator)->children[i]->M3);
#endif
            Tree_memory((*pbesurfaceoperator)->children[i]->M3);
         }

         unsigned int dielectriccavitychildcount = 0;

         for (j = 0; j < numdielectriccavities; j++) {
            if (dielectriccavityparent[j] != i)
               continue;

            printf("\nDielectric cavity surface %u is inside dielectric surface %u\n", j, i);
            printf("Constructing operator INSIDE dielectric cavity %u\n", j);

            SurfaceOperator dielectricparent = (*pbesurfaceoperator)->children[i];

            dielectricparent->children[dielectriccavitychildcount] = SurfaceOperator_allocate();
            dielectricparent->children[dielectriccavitychildcount]->kernel = POISSON_KERNEL;
            dielectricparent->children[dielectriccavitychildcount]->epsilon = outerdielectric;
            dielectricparent->children[dielectriccavitychildcount]->kappa = 0.0;
            dielectricparent->children[dielectriccavitychildcount]->mynumpanels = numdielectriccavitypanels[j];
            dielectricparent->children[dielectriccavitychildcount]->numchildren = 0;
            dielectricparent->children[dielectriccavitychildcount]->children = NULL;
            dielectricparent->children[dielectriccavitychildcount]->parent = dielectricparent;
            dielectricparent->children[dielectriccavitychildcount]->charges = NULL;

            dielectricparent->children[dielectriccavitychildcount]->resultInternal = Vector_allocate(numdielectriccavitypanels[j]);
            dielectricparent->children[dielectriccavitychildcount]->resultExternal = Vector_allocate(numdielectriccavitypanels[j]);

            Vector3D* dielectriccavitycentroids = (Vector3D*)calloc(numdielectriccavitypanels[j], sizeof(Vector3D));

            for (p = 0; p < numdielectriccavitypanels[j]; p++)
               dielectriccavitycentroids[p] = dielectriccavitypanels[j][p]->centroid;

            dielectricparent->children[dielectriccavitychildcount]->tree = Tree_allocate(dielectriccavitypanels[j], numdielectriccavitypanels[j], dielectriccavitycentroids, numdielectriccavitypanels[j], MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
            Tree_lists(dielectricparent->children[dielectriccavitychildcount]->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
            Tree_fill(dielectricparent->children[dielectriccavitychildcount]->tree);
#endif
            Tree_memory(dielectricparent->children[dielectriccavitychildcount]->tree);

            dielectriccavitychildcount++;
         }
      }
   }

   unsigned int index = 0;

   SurfaceOperator_initStartIndices((*pbesurfaceoperator), &index);

   if (index != (2*numtotalpanels)) {
      error("SurfaceOperator indexing is WRONG! (%u != %u)", index, 2*numtotalpanels);
   }
}

void generateQualocationOperator(QualocationOperator* qualocationoperator,
											PDBentry* PDBentries, unsigned int numPDBentries,
											Panel** dielectricpanels, unsigned int* numdielectricpanels, unsigned int numdielectrics,
											Panel** dielectriccavitypanels, unsigned int* numdielectriccavitypanels, unsigned int numdielectriccavities, unsigned int* dielectriccavityparent,
											unsigned int numtotalpanels) {
  Panel* allpanels = (Panel*)calloc(numtotalpanels, sizeof(Panel));
  Vector3D* allcentroids = (Vector3D*)calloc(numtotalpanels, sizeof(Vector3D));
  unsigned int i, j, count = 0;
  for (i = 0; i < numdielectrics; i++)
	 for (j = 0; j < numdielectricpanels[i]; j++) {
         allpanels[count] = dielectricpanels[i][j];
         allcentroids[count] = dielectricpanels[i][j]->centroid;
         count++;
	 }
  
  for (i = 0; i < numdielectriccavities; i++) {
	 if (dielectriccavityparent[i] == 999)
		continue;
	 
	 for (j = 0; j < numdielectriccavitypanels[i]; j++) {
		allpanels[count] = dielectriccavitypanels[i][j];
		allcentroids[count] = dielectriccavitypanels[i][j]->centroid;
		count++;
	 }
  }
  
  
  generateQualocationOperatorMinimal(qualocationoperator, PDBentries, numPDBentries, dielectricpanels, numdielectricpanels,
												 numdielectrics, dielectriccavitypanels, numdielectriccavitypanels, numdielectriccavities,
												 dielectriccavityparent, numtotalpanels);
  
  printf("\nConstructing qualocation operator for all panels\n");
  
  Tree_lists((*qualocationoperator)->tree);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
  Tree_fill((*qualocationoperator)->tree);
#endif
  Tree_memory((*qualocationoperator)->tree);
}

void generateQualocationOperatorMinimal(QualocationOperator* qualocationoperator,
PDBentry* PDBentries, unsigned int numPDBentries,
Panel** dielectricpanels, unsigned int* numdielectricpanels, unsigned int numdielectrics,
Panel** dielectriccavitypanels, unsigned int* numdielectriccavitypanels, unsigned int numdielectriccavities, unsigned int* dielectriccavityparent,
unsigned int numtotalpanels) {
   unsigned int a, i, j, count = 0;

   if (ionicstrength != 0.0)
      error("Qualocation method only works without salt!");

   *qualocationoperator = QualocationOperator_allocate();

   (*qualocationoperator)->innerdielectric = innerdielectric;
   (*qualocationoperator)->outerdielectric = outerdielectric;

   Panel* allpanels = (Panel*)calloc(numtotalpanels, sizeof(Panel));
   Vector3D* allcentroids = (Vector3D*)calloc(numtotalpanels, sizeof(Vector3D));

   for (i = 0; i < numdielectrics; i++)
      for (j = 0; j < numdielectricpanels[i]; j++) {
         allpanels[count] = dielectricpanels[i][j];
         allcentroids[count] = dielectricpanels[i][j]->centroid;
         count++;
      }

   for (i = 0; i < numdielectriccavities; i++) {
      if (dielectriccavityparent[i] == 999)
         continue;

      for (j = 0; j < numdielectriccavitypanels[i]; j++) {
         allpanels[count] = dielectriccavitypanels[i][j];
         allcentroids[count] = dielectriccavitypanels[i][j]->centroid;
         count++;
      }
   }
	printf("numPDBentries = %d\n", numPDBentries);
   (*qualocationoperator)->charges = Charge_allocate();
   (*qualocationoperator)->charges->numcharges = numPDBentries;
   (*qualocationoperator)->charges->points = (Vector3D*)calloc(numPDBentries, sizeof(Vector3D));
   (*qualocationoperator)->charges->charges = Vector_allocate(numPDBentries);

   for (a = 0; a < numPDBentries; a++) {
      Vector3D charge = Vector3D_allocate();
      charge->x = PDBentries[a].x;
      charge->y = PDBentries[a].y;
      charge->z = PDBentries[a].z;
      (*qualocationoperator)->charges->points[a] = charge;
      if (toupper(PDBentries[a].chain) != 'X')
         (*qualocationoperator)->charges->charges[a] = PDBentries[a].charge;
      else
         (*qualocationoperator)->charges->charges[a] = 0.0;
   }

  (*qualocationoperator)->tree = Tree_allocate(allpanels, numtotalpanels, allcentroids, numtotalpanels, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, DOUBLE_LAYER_INT, 0.0);

   if (accelerateM3) {
	  generateQualocationOperatorA1A3(*qualocationoperator);
   }
}

void generateQualocationOperatorA1A3(QualocationOperator qualocationoperator) {
  printf("\nConstructing combined M1 / M3 operator\n");
  unsigned int numtotalpanels = qualocationoperator->tree->numpanels;
  Panel *allpanels = qualocationoperator->tree->panels;
  qualocationoperator->M1M3 = Tree_allocate(allpanels, numtotalpanels, qualocationoperator->charges->points, qualocationoperator->charges->numcharges, MAX_PANELS_PER_FINEST_CUBE, POISSON_KERNEL, NULL, GRID_SIZE, SVD_ERROR, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
  Tree_lists(qualocationoperator->M1M3);
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
  Tree_fill(qualocationoperator->M1M3);
#endif
  Tree_memory(qualocationoperator->M1M3);
}

void generateSurfaceOperatorPreconditioner(Preconditioner* preconditioner, SurfaceOperator pbesurfaceoperator, unsigned int numtotalpanels) {
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
   *preconditioner =  Preconditioner_allocate(2 * numtotalpanels, 2 * numtotalpanels);
   SurfaceOperator_generatePreconditioner(pbesurfaceoperator, *preconditioner);
   //SurfaceOperator_generatePreconditioner_diagonal(pbesurfaceoperator, *preconditioner);
   //SurfaceOperator_generatePreconditioner_identity(pbesurfaceoperator, *preconditioner);
   Preconditioner_factor(*preconditioner);
#endif	
}

void generateQualocationOperatorPreconditioner(Preconditioner* preconditioner, QualocationOperator pbequalocationoperator, unsigned int numtotalpanels) {
#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
   *preconditioner =  Preconditioner_allocate(numtotalpanels, numtotalpanels);
   QualocationOperator_generatePreconditioner(pbequalocationoperator, *preconditioner);
   Preconditioner_factor(*preconditioner);
#endif
}

void generateRHS(Vector* rhs, SurfaceOperator pbesurfaceoperator, unsigned int numtotalpanels) {
   printf("There might be a memory leak in here...\n");
   *rhs = Vector_allocate(2 * numtotalpanels);

   SurfaceOperator_makeRHS(pbesurfaceoperator, *rhs);
}

void generateQualocationRHS(Vector* rhs, QualocationOperator pbequalocationoperator, unsigned int numtotalpanels) {
   *rhs = Vector_allocate(numtotalpanels);

   QualocationOperator_makeRHS(pbequalocationoperator, *rhs);   
}

void solveSurfaceOperator(Vector* solution, SurfaceOperator pbesurfaceoperator, Preconditioner preconditioner, Vector rhs, unsigned int numtotalpanels) {
   *solution = Vector_allocate(2 * numtotalpanels);

   GMRES_SurfaceOperator(pbesurfaceoperator, preconditioner, rhs, *solution, 2 * numtotalpanels, errortolerance);
}

void solveQualocationOperator(Vector* solution, QualocationOperator pbequalocationoperator, Preconditioner preconditioner, Vector rhs, unsigned int numtotalpanels) {
   *solution = Vector_allocate(numtotalpanels);

   GMRES_QualocationOperator(pbequalocationoperator, preconditioner, rhs, *solution, numtotalpanels, errortolerance);
}

void collectPotentials(PDBentry* PDBcollectentries, unsigned int numPDBcollectentries, PDBentry* PDBinputentries, unsigned int numPDBinputentries, SurfaceOperator pbesurfaceoperator, Vector solution) {
   Vector potentials = Vector_allocate(numPDBcollectentries);

   SurfaceOperator_calculateReactionPotentials(pbesurfaceoperator, solution, potentials);

   unsigned int a, c;

   for (a = 0; a < numPDBcollectentries; a++)
      PDBcollectentries[a].potential = KT_CONVERSION * potentials[a] / innerdielectric;

   if (addcoulombic) {
      for (c = 0; c < numPDBcollectentries; c++)
         for (a = 0; a < numPDBinputentries; a++)
            if (toupper(PDBinputentries[a].chain) != 'X') {
               real distance = sqrt((PDBinputentries[a].x - PDBcollectentries[c].x) *
                                    (PDBinputentries[a].x - PDBcollectentries[c].x) +
                                    (PDBinputentries[a].y - PDBcollectentries[c].y) *
                                    (PDBinputentries[a].y - PDBcollectentries[c].y) +
                                    (PDBinputentries[a].z - PDBcollectentries[c].z) *
                                    (PDBinputentries[a].z - PDBcollectentries[c].z));
               if (distance > 0.0)
                  PDBcollectentries[c].potential += KT_CONVERSION * PDBinputentries[a].charge / (innerdielectric * distance);
            }
   }

   Vector_free(potentials);
}

void collectQualocationPotentials(PDBentry* PDBcollectentries, unsigned int numPDBcollectentries, PDBentry* PDBinputentries, unsigned int numPDBinputentries, QualocationOperator pbequalocationoperator, Vector solution) {
   Vector potentials = Vector_allocate(numPDBcollectentries);

   QualocationOperator_calculateReactionPotentials(pbequalocationoperator, solution, potentials);

   unsigned int a, c;

   for (a = 0; a < numPDBcollectentries; a++)
      PDBcollectentries[a].potential = KT_CONVERSION * potentials[a] / innerdielectric;

   if (addcoulombic) {
      for (c = 0; c < numPDBcollectentries; c++)
         for (a = 0; a < numPDBinputentries; a++)
            if (toupper(PDBinputentries[a].chain) != 'X') {
               real distance = sqrt((PDBinputentries[a].x - PDBcollectentries[c].x) *
                                    (PDBinputentries[a].x - PDBcollectentries[c].x) +
                                    (PDBinputentries[a].y - PDBcollectentries[c].y) *
                                    (PDBinputentries[a].y - PDBcollectentries[c].y) +
                                    (PDBinputentries[a].z - PDBcollectentries[c].z) *
                                    (PDBinputentries[a].z - PDBcollectentries[c].z));
               if (distance > 0.0)
                  PDBcollectentries[c].potential += KT_CONVERSION * PDBinputentries[a].charge / (innerdielectric * distance);
            }
   }

   Vector_free(potentials);
}

void getPanelsFromArray(unsigned int numpanels,
								real *paneldata,
								Panel *panels)
{
  unsigned int i,j;
  Vector3D v1, v2, v3;
  v1 = Vector3D_allocate();
  v2 = Vector3D_allocate();
  v3 = Vector3D_allocate();
  for (i = 0; i < numpanels; i++) {
	 v1->x = paneldata[i*9+0]; v1->y = paneldata[i*9+1]; v1->z = paneldata[i*9+2];
	 v2->x = paneldata[i*9+3]; v2->y = paneldata[i*9+4]; v2->z = paneldata[i*9+5];
	 v3->x = paneldata[i*9+6]; v3->y = paneldata[i*9+7]; v3->z = paneldata[i*9+8];
	 FlatPanel fp = FlatPanel_allocate(v1, v2, v3);
	 panels[i] = Panel_allocate();
	 Panel_FlatPanel(panels[i], fp);
  }
  Vector3D_free(v1);
  Vector3D_free(v2);
  Vector3D_free(v3);
}

void solveBEMproblem(unsigned int *numpanels, real *paneldata,
							real *innerdielectric, real *outerdielectric,
							real *rhs, real *sol)
{
  unsigned int i,j;
  Panel *panels = (Panel *)calloc(*numpanels, sizeof(Panel));
  getPanelsFromArray(*numpanels, paneldata, panels);

  Vector3D *centroids = (Vector3D*)calloc(*numpanels, sizeof(Vector3D));
  for (i = 0; i < *numpanels; i++) {
	 centroids[i] = Vector3D_allocate();
	 Vector3D_copy(centroids[i], panels[i]->centroid);
  }

  real maxpanelsperfinestcube = 27;
  unsigned int gridpoints = 3;
  real svdtol = 1e-5;
  real gmrestol = 1e-5;
  
  Tree tree = Tree_allocate(panels, *numpanels, centroids, *numpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, DOUBLE_LAYER_INT, 0.0);
  Tree_lists(tree);
  Tree_fill(tree);

  Preconditioner preconditioner = Preconditioner_allocate(*numpanels, *numpanels);
  Preconditioner_fill_diagonal_solv_ecf_qual(preconditioner, tree,
															*innerdielectric, *outerdielectric);
  Preconditioner_factor(preconditioner);

  Vector Vec_rhs = Vector_allocate(*numpanels);
  Vector Vec_sol = Vector_allocate(*numpanels);
  for (i = 0; i < *numpanels; i++) { 
	 Vec_rhs[i] = rhs[i];
  }

  GMRES_solv_ecf_qual(tree, preconditioner, Vec_rhs, Vec_sol, gmrestol,
							 *innerdielectric, *outerdielectric);

  for (i = 0; i < *numpanels; i++) { 
	 sol[i] = Vec_sol[i];
  }

  // clean up, go home
  Vector_free(Vec_sol);
  Vector_free(Vec_rhs);
  Preconditioner_free(preconditioner);
  Tree_free(tree);
  for (i = 0; i < *numpanels; i++) {
	 Vector3D_free(centroids[i]);
	 Panel_free(panels[i]);
  }
  free(centroids);
  free(panels);
}
