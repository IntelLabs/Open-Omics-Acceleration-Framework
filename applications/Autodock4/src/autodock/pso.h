#ifndef _pso_h
#define _pso_h

#include "gs.h"
#include "ls.h" 
#include "structs.h"

extern FILE *logFile;

// define class for PSO algorithmic options
// to avoid having to change dozens of method signatures whenever
// the option set is changed. MP TSRI 2011
struct PSO_Options {
	double pso_w;	   // inertia weight
	double pso_w_start;	// pso_w at beginning of run
	double pso_w_end;	// pso_w at conclusion of run, see pso.cc
	double c1;	// cognitive
	double c2;	// social
	int pso_K;      // number of neighbor particles
	double c;    // constriction factor for cPSO   MP TODO notused
	double pso_vmax_scale; // MP 
	Boole pso_neighbors_dynamic; // MP
	Boole pso_neighbors_symmetric; // MP
	Boole pso_random_by_dimension; // MP
	Boole pso_adaptive_velocity; // MP
	Boole pso_regenerate_at_limit; // MP
	Boole pso_stage2constriction; // MP
	Boole pso_interpolate_as_scalars; // MP nothing else is implemented yet
  // default values for PSO options :
  public:
    PSO_Options () :
        pso_w(1.0), // w
        pso_w_start(0.9), // pso_w at beginning of run
        pso_w_end(0.4), // pso_w at conclusion of run
        c1(2.05), 
        c2(2.05),
        pso_K(4),
        c(0.01),
	pso_vmax_scale(0.1), // MP not yet implemented
        pso_neighbors_dynamic(false), // 
        pso_neighbors_symmetric(false), //  MP not yet implemented
        pso_random_by_dimension(true),  // 
        pso_adaptive_velocity(false),  //
        pso_regenerate_at_limit(true),  //
        pso_stage2constriction(false), //
	pso_interpolate_as_scalars(true)  //
        { }
	};

class ParticleSwarmGS : public Global_Search 
{
	private:
		Population *_Pi;	// best solution for each individual in its own searching history
		Individual	*_Pg;	// current best solution
		int best; // index in Pi of current global best solution
		unsigned int pop_size;	// population size
		int size;	// number of dimensions (7*nlig + num_torsion)
		float **v;	// velocity
		float *vmax;	//max velocity 
		float *vmin;	// min velocity
		double *xmax;	// max coord bound
		double *xmin;	// min coord bound
		PSO_Options pso_options;
	    
		unsigned int generations;
		int outputEveryNgens;
        Output_pop_stats output_pop_stats;
	
		Local_Search *LocalSearchMethod;
	
	public:	
		~ParticleSwarmGS();
		ParticleSwarmGS(
			float *init_vmax, 
			float *init_vmin, 
			double *init_xmax, 
			double *init_xmin, 
			PSO_Options init_pso_options, 
			Local_Search *init_LS, 
			unsigned int init_max_evals, 
			unsigned int init_max_generations, 
			Output_pop_stats output_pop_stats); 			
		
		Individual& getBest();	

      void initialize(const unsigned int, const unsigned int, const int, FILE *);
      unsigned int num_generations(void) const;
		
		// The following part are derived virtual functions
        char *shortname(void);
        char *longname(void);
		void reset(void);
        void reset(const Output_pop_stats&);
        int terminate(void);
        int search(Population &, Eval *, int outlev, FILE * logFile); 
	int localsearch(Population &, Local_Search *, Eval *evaluate, int outlev, FILE * logFile);
};

inline char * ParticleSwarmGS::shortname(void)
{
        return "PSO";
}

inline char * ParticleSwarmGS::longname(void)
{
        return "PARTICLE SWARM OPTIMIZATION";
}

inline ParticleSwarmGS::~ParticleSwarmGS()
{
	if(_Pi)
		delete _Pi;
	if(_Pg)
		delete _Pg;
	if(v) {
		for(unsigned int i=0;i < pop_size; i++)	
			delete [] v[i];
		delete [] v;
	}
}

inline ParticleSwarmGS::ParticleSwarmGS(
			float *init_vmax, 
			float *init_vmin, 
			double *init_xmax, 
			double *init_xmin, 
			PSO_Options init_pso_options, 
			Local_Search *init_LS, 
			const unsigned int init_max_evals, 
			const unsigned int init_max_generations, 
			Output_pop_stats init_output_pop_stats) :
    Global_Search(  init_max_evals, init_max_generations)
{
	vmax = init_vmax; 
	vmin = init_vmin;			  
	xmax = init_xmax; 
	xmin = init_xmin;			  
	pso_options = init_pso_options;
	LocalSearchMethod = init_LS;
    generations = 0;
	output_pop_stats = init_output_pop_stats;
	_Pi =NULL; 
	_Pg = NULL ;
	 v = NULL; 
}

inline Individual& ParticleSwarmGS:: getBest()
{
	return *_Pg;
}
inline void ParticleSwarmGS::initialize(const unsigned int init_pop_size, const unsigned int ndims, int outlev, FILE *logFile)
{

    pop_size = init_pop_size;
    size = ndims;
}

// MP adapted from gs.h (can we just have one of these?)
inline unsigned int ParticleSwarmGS::num_generations(void) const
{
   return(generations);
}

// The following part are derived virtual functions
//int search(Population &, Eval *evaluate, int outlev, FILE * logFile);

inline int ParticleSwarmGS::terminate(void)
{
   if (max_generations>(unsigned) 0) {
      return((unsigned)generations>=max_generations); 
   } else {
      return(0);  //  Don't terminate
   }
}

	
inline void ParticleSwarmGS::reset(void)
{
	generations = 0;
	//MP pso_w = pso_w_start;
	if(_Pi)
		delete _Pi;
	if(_Pg)
		delete _Pg;
	if(v) {
		for(unsigned int i=0;i < pop_size; i++)	
			delete [] v[i];
		delete [] v;
	}
	
	_Pi = NULL;
	_Pg = NULL;
	v = NULL;
}
	
inline void ParticleSwarmGS::reset(const Output_pop_stats &init_output_pop_stats)
{
    output_pop_stats = init_output_pop_stats; 
	reset();
}
	
/*
Constants used for the PSO Work
 */
#define PSO_D_MAX     (7+MAX_TORS)   // Max number of dimensions of the search space
#define PSO_S_MAX 	  200	 // Max swarm size
#define PSO_K_MAX 	  PSO_S_MAX	 // Max neighbors


/*______________________________________________________________________________
**PSO Work Structures */

typedef struct velocity {
	int size;
	double v[PSO_D_MAX];
} Velocity;

/*______________________________________________________________________________*/

typedef struct position {
        int size;
        double x[PSO_D_MAX];
        double f;		// fitness value of particle
        double prev_x[PSO_D_MAX];		// previous pos of particle
} Position;



#endif
