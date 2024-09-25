#include "pso.h"
#include "ranlib.h"
#include "support.h" // for evalmode and rep_constants
#include <math.h>

//extern Eval evaluate;
 // MPique - the following static vars should be put in the class...
static Boole init_links;
static float prevBestE;
static	float prevE[PSO_S_MAX];       // E array for particles in previous 
static	float curE[PSO_S_MAX];        // E array for particles in current step


inline float Norm(float *x, int n)
{
	float s = 0;
	for(int i=0;i<n;i++)
		s += (x[i]*x[i]);
	return sqrt(s);
}


/***********************************************************************
 * Particle Swarm Optimization (PSO) with time varying inertia weight w.
 * Global PSO and local S&W search.
 * 
 * Based on PSO in SODOCK: J. Comput Chem, 28, 2, 612-623, 2007.
 * 
 * Huameng Li, 07/2008
 * 
 *
 * Key modifications spring-summer 2011 by R Huey & M Pique at TRSI:
 *   No longer does local search here, the caller will do that
 *     after each generation using our localsearch method, below.
 * Caution: do not reorder the Pop array outside of these functions. MP/RH
 * Packaged algorithm options into a structure to simplify signatures.
 *
 * Caution: not yet thread-safe (local Pi array) and does not delete
 *  the "new"-created arrays, which are members of the ParticleSwarmGS class, see pso.h - M Pique 2014
 ***********************************************************************/
int ParticleSwarmGS::search(Population &Pop, Eval *evaluate, int outlev, FILE * logFile)
{
	int j;
	double curVal;
	double newVal;
	double Pop_best_value;
	float g_ratio, e_ratio, ratio;
	int links[PSO_S_MAX][PSO_K_MAX]; // neighbors of i (the [i][0] is always self)


	// copy some options for convenience
	float pso_w = pso_options.pso_w;	   // inertia weight
	float w_start = pso_options.pso_w_start;	// pso_w at beginning of run
	float w_end = pso_options.pso_w_end;	// pso_w at conclusion of run
	float c1 = pso_options.c1;
	float c2 = pso_options.c2;
	int pso_K = pso_options.pso_K;      // number of neighbor particles

	
	// on first call per run, allocate Pi array, allocate and initialize velocity vectors, prevE, curE
	// print pso options if outlev>0
	if(_Pi == NULL) {
#define bool(i)  ((i)?"true":"false")
		
		pr(logFile, "PSO max_generations = %d\n", max_generations);
		pr(logFile, "PSO max_evaluations = %d\n", max_evals);
		pr(logFile, "PSO size = %d\n", size);

		if(outlev>LOGBASIC){
		pr(logFile, "PSO w_start = %.2f\n", pso_options.pso_w_start);
		pr(logFile, "PSO w_end = %.2f\n", pso_options.pso_w_end);
		pr(logFile, "PSO c1 = %.2f\n", pso_options.c1);
		pr(logFile, "PSO c2 = %.2f\n", pso_options.c2);
		pr(logFile, "PSO K = %d\n", pso_options.pso_K);
		pr(logFile, "PSO neighbors_dynamic = %s\n", bool(pso_options.pso_neighbors_dynamic));
		pr(logFile, "PSO neighbors_symmetric = %s\n", bool(pso_options.pso_neighbors_symmetric));
		pr(logFile, "PSO random_by_dimension = %s\n", bool(pso_options.pso_random_by_dimension));
		pr(logFile, "PSO adaptive_velocity = %s\n", bool(pso_options.pso_adaptive_velocity));
		pr(logFile, "PSO regenerate_at_limit = %s\n", bool(pso_options.pso_regenerate_at_limit));
		pr(logFile, "PSO stage2constriction = %s\n", bool(pso_options.pso_stage2constriction));
		pr(logFile, "PSO interpolate_as_scalars = %s\n", bool(pso_options.pso_interpolate_as_scalars));
		}
				
		fprintf(logFile, "PSO Allocate initial velocity of particles...\n");
		_Pi = new Population(Pop); // copy Pop
		
		pop_size = Pop.num_individuals();
		pr(logFile, "PSO pop_size = %d\n", pop_size); // MP debug
						
		best = 0;
		for(unsigned int i = 1; i < pop_size; i++)
			if( (*_Pi)[i].value(Normal_Eval) < (*_Pi)[best].value(Normal_Eval) )
				best = i;
		_Pg = new Individual( (*_Pi)[best] );
		
		// allocate velocity 
		v = new float * [pop_size];
		for(unsigned int i=0; i < pop_size; i++) {
			v[i] = new float [size];
		}
			
		// initial velocity of particles
		for(unsigned int i = 0; i < pop_size; i++) {			
			for(j = 0; j < size; j++) {
				v[i][j] = random_range(vmin[j], vmax[j]);				
			}			
		}
		// MP: note that with adaptive velocity, the BIG will prevent the
		//  first velocity update from occurring
		for(unsigned int i = 0; i < pop_size; i++) prevE[i] = curE[i] = BIG; // initially unfavorable

		if(outlev>LOGRUNV) {
		pr(logFile, "PSO Initial velocity V:\n");
		for(unsigned int i=0;i<pop_size;i++) {
			pr(logFile, " V[%d] ", i);
			for(j=0;j<size;j++)
				pr(logFile, "%+6.2f ", v[i][j]);
			pr(logFile, "\n");
		}				
		pr(logFile, "Pi:\n");
		for(unsigned int i=0;i<pop_size;i++) {
			pr(logFile, "Pi[%d] ", i);
			for(j=0;j<size;j++)
				pr(logFile, "%+6.2f ", (*_Pi)[i].phenotyp.gread(j).real);
			pr(logFile, "\n");
		}
		}
								
		// Display field title
		//pr(logFile, "Generation NumEvals     Pg       PopAvg     PiAvg       w      |V| Avg\n");
		//pr(logFile, "---------- -------- ---------- ---------- ---------- -------- --------\n");					
		//pr(logFile, "Generation\tNumEvals\tpi Pop_best\tgBest Pg\tImproved\n");	
		pr(logFile, "Generation NumEvals  gBest E   \n");
		pr(logFile, "---------- -------- ---------- \n");					
		init_links = TRUE;
	}
	
	// Update weights as the run progresses
	// Huameng had:  ratio = (float)evaluate.evals() / num_evals;
	// TSRI (MP) uses larger of fraction of generations used & evaluations used...
	g_ratio = max_generations>0 ? (float)generations / max_generations : 0;
	e_ratio = max_evals>0 ? (float)evaluate->evals() / max_evals : 0;
	ratio =  max(g_ratio, e_ratio);
	if(pso_options.pso_stage2constriction && ratio>0.8) {
		// Huameng: Incorporate stage 2 constriction PSO 
		pso_w = 0.1;
		c1 = 6.05;
		c2 = 6.05;
		}
	else pso_w = w_end +(w_start - w_end) * (1.0 - ratio);
	
	// Set shorthand names for: best for each indiv, best individual in pop, best in neighborhood
	Population &Pi = (Population &)(*_Pi);
	Individual &Pg = (Individual &)(*_Pg);
	Individual &nbBest = (Individual &)(*_Pg);
	
	//double piBestE[pop_size];	    // E array for personal Best	

	// set prevE from curE (previous move)
	for(unsigned int i=0;i<pop_size;i++) prevE[i] = curE[i];

	// if global best has not improved, recreate neighborhood links before this move
	if(pso_options.pso_neighbors_dynamic && Pg.value(Normal_Eval) > prevBestE) {
		if(outlev>LOGRUNV)  {
			// DEBUG MP
			fprintf(logFile, "PSO init_links Pg= %.2f prevBestE= %.2f\n",
			 Pg.value(Normal_Eval) , prevBestE);
			}
		init_links=TRUE;
	}

	for(unsigned int i = 0; i < pop_size; i++ ) {		
		curE[i] = Pop[i].value(Normal_Eval);
		// update Pi from Pop
		if( Pop[i].value(Normal_Eval) < Pi[i].value(Normal_Eval) )		
		        Pi[i] = Pop[i];
		 // assert energy of Pi[i] <= energy of Pop[i]
	}
	// set Global Best Pg, which points to _Pg
	best = 0;	
	for(unsigned int i = 1; i < pop_size; i++ ) {
		if( Pi[i].value(Normal_Eval) < Pi[best].value(Normal_Eval) ) {
			best = (int)i;			
		}
	}
	Pg = Pi[best];
	prevBestE=Pg.value(Normal_Eval);
	// assert energy of Pg <= energy of Pi[*] <= energy of Pop[*]
		   		   	 	     	   	 	   	   	   	  
	// set up or modify neighborhood lists
	//  dynamic option remakes neighborhoods if Pi not improved
	//     this dynamic neighborhood is from the "other" PSO code
	// if the entire population is the neighborhood, this table is not used
        if (init_links && pso_K<(int)pop_size) {
            //Who informs who, at random
	    // MP:note these are asymmetric links - unlike Huamengs I think
            for (unsigned int s=0; s<pop_size; s++) {
                links[s][0]=s; // Each particle informs itself
	        // K-1 other links (could possibly be self or repeated)
                for(int i=1; i<pso_K;i++) links[s][i]= (int) random_range(0, pop_size-1);
            }        
	    init_links=FALSE;
        }


	// compute velocity vector for each particle
	// Note - this loop uses Pop[i] but does not modify Pop[i]
	for(unsigned int i = 0; i < pop_size; i++) {
				
		if(pso_K < (int) pop_size) {
			// get the best in neighborhood of "i", e.g., Neighborhood Best (nbBest)		
			int g = i; // self
			for(j =  1; j < pso_K; j++) {
				if(Pi[links[i][j]].value(Normal_Eval) < Pi[g].value(Normal_Eval))
				g = links[i][j];
			}
			nbBest=Pi[g];
		}
	        else nbBest=Pg;  // use global best as neighborhood best

				
		double r1, r2;
		r1 = ranf();
		r2 = ranf();

		// size is the dimension of docking search (e.g, n*nlig + ntor)			
		for(j = 0; j < size; j++) {
			if(pso_options.pso_random_by_dimension) {
				r1 = ranf();
				r2 = ranf();
			}

			// note that phenotype "gread" is x,y,z,qx,qy,qz,qw,t1...
			curVal = Pop[i].phenotyp.gread(j).real;
									
			
			// update velocity

			// if not adaptive_velocity, always update velocity.
			// if adaptive_velocity, update velocity only if individual's energy
			//  is WORSE than (or same as) its previous energy
			if( (!pso_options.pso_adaptive_velocity) || Pop[i].value(Normal_Eval)>=prevE[i]) {
				v[i][j] = pso_w * v[i][j] 
					+ c1 * r1 * (Pi[i].phenotyp.gread(j).real - curVal)
					+ c2 * r2 * (nbBest.phenotyp.gread(j).real - curVal);
				}
													   // MP TODO handle quat and torsion velocities better
			                      																			
			// verify and restrict the new velocity component
			if(v[i][j] > vmax[j]) v[i][j] = vmax[j];
			else if(v[i][j] < vmin[j]) v[i][j] = vmin[j];
		} // next dim
	}	// end compute velocity 
				
	// update position (translation, rotation, torsions) Pop[i] of each particle			
	for(unsigned int i = 0; i < pop_size; i++) {
		for(j = 0; j < size; j++) {

			curVal = Pop[i].phenotyp.gread(j).real;
			// wrap torsions (MP not sure if ever necessary...)
			if(is_conformation_index(j)) curVal=WrpModRad(curVal);

			// apply velocity
			newVal = curVal + v[i][j];


			// wrap torsions
			// and verify new translation value is within bounds
			// MP: note that only translation has a bound -
			//  quaternions and torsions wrap 
			//  (for quaternions, in ways we are ignoring here).
			if(is_conformation_index(j)) newVal=WrpModRad(newVal);
			else if( is_translation_index(j)
			    && (newVal < xmin[j] || newVal > xmax[j]) ) {
				if(pso_options.pso_regenerate_at_limit) {
				// regenerate at random xyz - do not 
				//  change orientation or torsions for now
				if(outlev>LOGRUNV)  fprintf(logFile, 
				"PSO [%d] hit limit[%c] (%.2f,%.2f,%.2f) ", i, "XYZ"[j],
				Pop[i].phenotyp.gread(0).real,
				Pop[i].phenotyp.gread(1).real,
				Pop[i].phenotyp.gread(2).real);

				for(int d=0;d<3;d++) {
					Pop[i].phenotyp.write( random_range(xmin[d], xmax[d]), d);
					}
				if(outlev>LOGRUNV)  fprintf(logFile, "regen at (%.2f,%.2f,%.2f)\n",
				Pop[i].phenotyp.gread(0).real,
				Pop[i].phenotyp.gread(1).real,
				Pop[i].phenotyp.gread(2).real);
				break;  // end of this individual move
				}
			else {
				// if not regenerate at random,
				// simply clamp this translation component
			     if(newVal < xmin[j]) newVal=xmin[j];
			     else if (newVal > xmax[j]) newVal=xmax[j];
				}
			}
													// MP TODO handle quat and torsion accumulation better
										
			// update x,y,z, quaternion, torsion of particle i
			Pop[i].phenotyp.write(newVal, j);
		} // next dim
	// must always normalize Quaternion after modifying its components
	Quat q = Pop[i].phenotyp.readQuat();
	Pop[i].phenotyp.writeQuat( normQuat( q ));

	Pop[i].inverse_mapping(); // MP@@ copy phenotype to genotype (may not be necc)
		
	}	// end update position
	
	// Update personal best Pi after this swarm search generation.
	// Find the best in Pop after this swarm search generation.
	// Note each could be better or worse than former best, stored in Pg
	// This 'best' index will be used by a subsequent call to LocalSearch,
	//   see below. M Pique June 2011
	best = 0;
	Pop_best_value = BIG;
	for(unsigned int i = 0; i < pop_size; i++) {		
		double piCurE = Pop[i].value(Normal_Eval);
		if(piCurE < Pop_best_value) {
			best = i;
			Pop_best_value = piCurE;
		}
		
		// Update Pi, personal Best in history				
		if(piCurE < Pi[i].value(Normal_Eval)) {
			Pi[i] = Pop[i];			
		}										
	}

	generations++;
	
#ifdef SWARM_ANALYSIS
	// comments out the swarm analysis   -Huameng
	float Pop_avg = 0;
	float Pi_avg = 0;
	float v_avg = 0;
	for(i = 0; i < pop_size;i++) {
		Pop_avg += Pop[i].value(Normal_Eval);
		Pi_avg += Pi[i].value(Normal_Eval);
		v_avg += Norm(v[i], size);
	}
	Pop_avg /= pop_size;
	Pi_avg /= pop_size;
	v_avg /= pop_size;
		
	float diff[50];
	float dist = 0;
	int k;	
	int count = 0;
	for(i = 0; i < pop_size-1; i++) {
		for(j = i+1;j < pop_size; j++)
			for(k = 0; k < size; k++)
				diff[k] = Pop[i].phenotyp.gread(k).real - Pop[j].phenotyp.gread(k).real;
		dist += Norm(diff, size);
		count++;
	}
	dist /= count;
#endif /* SWARM_ANALYSIS */
	
	// Output PSO statistics
	if(outputEveryNgens > 0 && 
	  (generations % outputEveryNgens == 0||generations==1)) {
		//pr(logFile, "%d %8d %10.2f %10.2f %10.2f %6.2f %8.2f\n", generations, evaluate->evals(), Pg.value(Normal_Eval), Pop_avg, Pi_avg, w, v_avg);
		//ORIG pr(logFile, "%8d %10ld %10.2f  \n", generations, evaluate->evals(), Pg.value(Normal_Eval));		
		fflush(logFile);
	}
	//pr(logFile, "%8d\t%6d\t%6.3f\t%6.3f\t%6.3f\n", generations, evaluate->evals(), Pop_best_value, Pg.value(Normal_Eval), dist);	
	//pr(logFile, "%6d\t%6d\t%6.3f\t%6.3f\n\n", generations, evaluate->evals(), Pop_best_value, Pg.value(Normal_Eval));		
	//fflush(logFile);
	
	return (0);
}

int ParticleSwarmGS::localsearch(Population &Pop, Local_Search *local_method, Eval *evaluate, int outlev, FILE * logFile)
{
// does local search (unconditonally) on "best" (index) indiv, which may not be historic best

// Set shorthand names for best for each indiv, best individual in pop
Population &Pi = (Population &)(*_Pi);
Individual &Pg = (Individual &)(*_Pg);
	if(local_method == NULL) return(-1);
	if(best<0||(unsigned)best>=Pop.num_individuals()) {
		stop("PSO LocalSearch bad best value");
	}
	// MP DEBUG 2011
	if(outlev>LOGRUNV) (void)fprintf( logFile, "PSO LS only Pop[best=%d] was %f\n", 
	 best, Pop[best].value(Normal_Eval));

	// If local seach method is defined, apply local search			
	local_method->search(Pop[best], evaluate, outlev, logFile);

	if(Pop[best].value(Normal_Eval) < Pi[best].value(Normal_Eval)) {
		Pi[best] = Pop[best];			
		Pg = Pi[best]; // mp maybe not needed
	}					


	// MP DEBUG 2011:
	if(outlev>LOGRUNV) (void)fprintf( logFile, "PSO LS only Pop[best=%d] now %f\n", 
	 best, Pop[best].value(Normal_Eval));
   if (outlev > LOGRUNV) {
	(void)fprintf( logFile, " %f",Pop[best].value(Normal_Eval)); 
	(void)fprintf( logFile, " \n"); 
	    }
	return(0);
}
