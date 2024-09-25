/* main.cc */
#include <stdio.h>
#include <math.h>
#include <time.h>      // time_t time(time_t *tloc);
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <ctype.h> // tolower
#include <vector>

#include "call_cpso.h"
#include "printdate.h"
#include "timesyshms.h"
#include "assert.h"
#include "eval.h"
#include "qmultiply.h"
#include "ranlib.h"
#include "support.h"
#include "ls.h"


extern FILE *logFile;
extern Eval evaluate;

Representation **cnv_state_to_rep2(const State &state)
{
   register int i;
   Representation **retval;

   retval = new Representation *[5];
   retval[0] = new RealVector(1);
   retval[0]->write(state.T.x, 0);
   retval[1] = new RealVector(1);
   retval[1]->write(state.T.y, 0);
   retval[2] = new RealVector(1);
   retval[2]->write(state.T.z, 0);
   retval[3] = new RealVector(4);
   retval[3]->write(state.Q.x, 0);
   retval[3]->write(state.Q.y, 1);
   retval[3]->write(state.Q.z, 2);
   retval[3]->write(state.Q.w, 3);
   retval[4] = new RealVector(state.ntor);
   for(i=0; i<state.ntor; i++) {
      retval[4]->write(state.tor[i], i);
   }

   return(retval);
}



Individual cnv_state_to_ind2(const State &original, const int ntor)
{
   // BEGIN ADDITION
   // Added by gmm, 27-MAR-97, to solve these compiler warnings:
   //
   // call_ls.cc:59: warning: In this statement, the initialization of a non-const reference requires a temporary for "Genotype(5,cnv_state_to_rep(original))". (reftemporary)
   // call_ls.cc:59: warning: In this statement, the initialization of a non-const reference requires a temporary for "Phenotype(5,cnv_state_to_rep(original))". (reftemporary)
   // call_ls.cc:59: warning: In this statement, the initialization of a non-const reference requires a temporary for "(Individual(Genotype(5,cnv_state_to_rep(original)),Phenotype(5,cnv_state_to_rep(original))))". (reftemporary)

   Genotype temp_Gtype;
   Phenotype temp_Ptype;

   temp_Gtype = Genotype(5, cnv_state_to_rep2(original));
   temp_Ptype = Phenotype(5, cnv_state_to_rep2(original));

   Individual temp(temp_Gtype, temp_Ptype);

   return(temp);
   // END ADDITION
}


State call_cpso(Local_Search *const local_method, 
                const State& sInit, 
                const int n_exec, const int S, //swarm size
                const int D, //number of dimensions:7 + ntor
                double *const xmin, double *const xmax, 
                const unsigned int num_evals, 
                const int K, ConstDouble c1, ConstDouble c2,
                const int outlev)
{
    int nb_eval = 0;
    Position Xi[S_max]; //Max swarm size set to 1024 in constants.h
    Velocity Vi[S_max];
    Position P[S_max];
    State sNew[S_max];
    State sTemp;
    Individual iTemp;
    int g, s, d, i, m , best;
    double pso_energy;
    double energy_prev;
    double min= 0.0;
    int init_links;
    int LINKS[S_max][S_max];
    float prev_value[S_max];
    int evaluations=0;
    double Vmin[D_max], Vmax[D_max];
    //Constants for the Constriction PSO 
    double chi,phi;
    //State sbNew[S_max];
    std::vector<State> sbNew(S_max, State(sInit));
    int sb = 0;
    int ntor = sInit.ntor;
    evaluate.reset();


    // Calculating chi for CPSO
    phi = c1 + c2;

    //Use formula from Clerc, M. Kennedy, J., Evolutionary
    //Computation, IEEE Transactions on , vol. 6, no. 1, 58-73 (2002).
    //TODO: check reference
    chi = 2.0 / (2.0 - phi - sqrt( (phi * phi) - (4 * phi) ));
    chi = fabs(chi);
    fprintf(logFile, "\nConstants used for the Constriction PSO \nc1: %lf c2: %lf phi: %lf chi: %lf S: %d\n", c1, c2, phi, chi, S);
    // Initialization of the Particles
    for( s = 0; s < S ; s++)
    {
        //Updating evaluations        
        evaluations++;
        if (outlev > LOGRUNBASIC )
        {
            pr(logFile, "PSO: Updating the swarm at move %d (= %d and %u evaluations)\n", nb_eval+1, evaluations, (unsigned int) evaluate.evals());
        }
        Xi[s].size = D; 
        Vi[s].size = D;

        //Initialise the Particle 
        initialiseParticle(s, D, Xi, Vi, xmin, xmax, Vmin, Vmax);
        sNew[s]=sInit;
        copyDimension(&sNew[s], Xi[s]);
                
        //First Evaluations
        Xi[s].f = evaluate.evalpso(&sNew[s]);//E-test
        if (outlev >= 2) {
            pr(logFile, "\n(Initialization) Particle %d has Energy= %8.2lf\n", s, Xi[s].f);
            printState(logFile, sNew[s], 0);
        }
        P[s] = Xi[s];// Best Postion = current One
        // RG store the current energy value of particle 
        prev_value[s] = BIG;
    } //for loop on index 
    //Find the best
    best=0;
    for (s=0; s<S; s++)
    {
        if (P[s].f < P[best].f) best=s;
    }

    pso_energy = P[best].f;
    if (n_exec==1) min=pso_energy;
    energy_prev=pso_energy;
    init_links = 1;
    nb_eval=0;
    //loop: ITERATIONS
    //Iteration of the runs
    //for (nb_eval=1; nb_eval < num_evals; nb_eval++)
    do
    {
        if (init_links == 1)
        {
            //Who informs who, at random
            for (s=0; s<S; s++)
            {
                for(m=0; m<S; m++) LINKS[m][s]=0; //Init to "no link"
                LINKS[s][s]=1; // Each particle informs itself
            }
            for (m=0; m<S; m++) //Other Links
            {
                for(i = 0; i < K; i++)
                {
                    s = (int) random_range(0, S-1);
                    LINKS[m][s]=1;
                }
            }        
        }
        sb = 0;
        // The Swarm moves
        for (s=0; s<S; s++)
        {
            //Updating evaluations
            evaluations++;
            if (outlev >= LOGRUNV)
            {
                pr(logFile, "\nPSO: Updating the swarm at move %d (= %d and %u evaluations)\n", nb_eval+1, evaluations, (unsigned int) evaluate.evals());
            }
            //..find the best informant
            g = s;
            for (m=0; m<S; m++)
            {
                if (LINKS[m][s] == 1 && P[m].f < P[g].f) g = m;
            }
            //..compute the new velocity and move
            for (d=0; d<D; d++)
            {
                // -----Constriction -PSO ------Begin-
                if (Xi[s].f >= prev_value[s]) // Relaxation Velocity Update
                {
                    Vi[s].v[d] = chi * (Vi[s].v[d] + c1 * random_range(0,1) * (P[s].x[d] - Xi[s].x[d]) + c2 * random_range(0,1) * (P[g].x[d] - Xi[s].x[d]));
                }
                if (Vi[s].v[d] > Vmax[d]) 
                {        
                    Vi[s].v[d] = Vmax[d]; 
                }

                Xi[s].x[d] = Xi[s].x[d] + Vi[s].v[d];
                // -----Constriction -PSO ------End---
            }
            //..interval confinement (keep in the box) 
            for (d=0; d<D; d++)
            {
                if (d>5)
                {
                    if (Xi[s].x[d] > PI)
                    {
                        //printf("\n Value Initial = %lf", Xi[s].x[d]);
                        Xi[s].x[d] = Xi[s].x[d] - ((int)(Xi[s].x[d] / PI) * PI);
                        //Vi[s].v[d]= 0.0 ; 
                    }
                    if (Xi[s].x[d] < -PI)
                    {
                        //Xi[s].x[d] %= -PI;
                        Xi[s].x[d] = Xi[s].x[d] - ((int)(Xi[s].x[d] / -PI) * -PI);
                        //Vi[s].v[d]= 0.0 ; 
                    }
                }
                else
                {
                    if(Xi[s].x[d] < xmin[d]) 
                    {
                        Xi[s].x[d]= random_range(xmin[d], xmax[d]); 
                        Vi[s].v[d]= random_range(Vmin[d], Vmax[d]) ; 
                    }
                    if(Xi[s].x[d] > xmax[d])
                    {
                        Xi[s].x[d]= random_range(xmin[d], xmax[d]); 
                        Vi[s].v[d]= random_range(Vmin[d], Vmax[d]) ; 
                    }
                }
            }//for
            // copy the dimension to State Variable 
            copyDimension(&sNew[s], Xi[s]);
            //...evaluate the new position
            prev_value[s] = Xi[s].f;
            Xi[s].f= evaluate.evalpso(&sNew[s]);//E-test
            if (outlev >= LOGRUNVV) 
            {
                pr(logFile,"\nSwarmMove: (%d) \tParticle:  %d \tEnergy= %8.2lf\n", nb_eval + 1, s + 1, Xi[s].f);
                printState(logFile, sNew[s], 0);
            }
            if(Xi[s].f <= Xi[sb].f) sb=s;
        } //s
        copyDimension(&sbNew[s], Xi[sb]);
        if (local_method !=NULL)
        {
            iTemp = cnv_state_to_ind2(sbNew[s], (D-7));
            local_method->search(iTemp);
            sbNew[s] = iTemp.state(ntor);
            //pr(logFile, "\nsbNew[%d].T.x==%f @@@\n", s, sbNew[s].T.x);
        };
        for (s=0; s<S; s++)
        {
            //... Update the best previous position
            if (Xi[s].f < P[s].f)
            {
                P[s] = Xi[s];
                //...update the best of the bests
                if(P[s].f < P[best].f) best = s;
            }
        }        
        nb_eval++;
        //Swarm - Begin Activity moved as a function 
        swarmActivity(S, D, Xi, nb_eval, outlev);
        //Swarm  - End
        // If no improvement, information links will be reinitialised
        pso_energy=P[best].f;
        if(pso_energy >= energy_prev) init_links = 1;
        else init_links = 0;
        energy_prev = pso_energy;
        if (outlev >= LOGRUNV) 
        {
            pr(logFile, "PSO- Run: %2d \tPSObest Energy@Swarm_Move: %4d \tP= %8.2lf  \t(nbeval= %d and %u )\n", n_exec+1, nb_eval+1, pso_energy, evaluations, (unsigned int) evaluate.evals());
        };
    } while(evaluate.evals() < num_evals);
    //Result is stored in the History state for writePDBQ and for clustering ....
    sTemp=sInit;
    copyDimension(&sTemp, P[best]);
    mkUnitQuat(&sTemp.Q); // MP TODO probably unnecc
    //printState(logFile, sTemp, 2);        
    return (sTemp);
} //call_cpso

/* END OF PROGRAM */
/* EOF */
