#ifndef COPYDIMENSION
#define COPYDIMENSION

//#include "main.h"
#include "structs.h"
#include "constants.h"
#include "alea.h"

//void initialiseState( State *S );
//void copyState( State *destination, State  source);

void copyDimension( /* not const */ State *const S, const Position& R);
void copyState2Dimension( /* not const */ Position *const R, const State& S);

void initialiseDimension(const GridMapSetInfo *const info, /* not const */ double *const xmin, /* not const */ double *const xmax, const int D);
//void initialiseDimension(float xlo, float xhi, float ylo, float yhi, float zlo, float zhi, double *xmin, double *xmax, int D);
void initialiseParticle(const int s, const int D, /* not const */ Position *const Xi, /* not const */ Velocity *const Vi, const double *const xmin, const double *const xmax, double *const Vmin, double *const Vmax);
void swarmActivity(const int S, const int D, const Position *const Xi, const int nb_eval, const int outlev);

/*void printState( FILE *fp,
		 State state, 
		 int detail );

void writeState( FILE *fp, 
		 State state );

int checkState(State *D);

Molecule copyStateToMolecule(State *source, Molecule *mol);*/

#endif
