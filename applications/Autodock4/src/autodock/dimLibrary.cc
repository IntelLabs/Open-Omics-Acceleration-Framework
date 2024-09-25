/* dimLibrary.cc */

#include <math.h>

#include <stdio.h>
#include "pso.h"
#include "dimLibrary.h"
#include "constants.h"
#include "ranlib.h"

extern FILE *logFile;

// Modified May 2011 MP TSRI to use quaternion fields rather than axis-angle

void copyDimension( /* not const */ State *const S, const Position& R)
{
        register int i, j;
        S->T.x = R.x[0];
        S->T.y = R.x[1];
        S->T.z = R.x[2];
        S->Q.x = R.x[3];
        S->Q.y = R.x[4];
        S->Q.z = R.x[5];
        S->Q.w = R.x[6];
        for(j=7, i=0; i<S->ntor; i++, j++)
        {
		S->tor[i] = R.x[j];
	}
}


void copyState2Dimension(Position *const R , const State& S)
{
        register int i, j;
	R->x[0] = S.T.x;
	R->x[1] = S.T.y;
        R->x[2] = S.T.z;
        R->x[3] = S.Q.x;
        R->x[4] = S.Q.y;
        R->x[5] = S.Q.z;
        R->x[6] = S.Q.w;
        for(j=7, i=0; i<S.ntor; i++, j++)
        {
		R->x[j] = S.tor[i];
        }
}


void initialiseDimension(const GridMapSetInfo *const info,  /* not const */ double *const xmin, /* not const */ double *const xmax, const int D)
{
	int d;

    for (d=0; d < D ; d++)
        {
            
        	switch(d)
            {
            case 0: case 1: case 2:
				xmin[d] = info->lo[d];
                xmax[d] = info->hi[d];
				break;

			case 3: 
			case 4: 
			case 5: 
			case 6: 
			// For Quaternion terms
				xmin[d] = -1;
				xmax[d] = 1;
				break;
				
			default: // For Torsional Angles
				xmin[d] = -PI;
				xmax[d] = PI;
				break;
				
		}
	}
}


void initialiseParticle(const int s, const int D, /* not const */ Position *const Xi, /* not const */ Velocity *const Vi, const double *const xmin, const double *const xmax, double *const Vmin, double *const Vmax)
{
        int d;
	double temp;
    //printf("in initialiseParticle: s=%d, D=%d, Xi=%f  \n",s,D,*Xi);
    //printf("in initialiseParticle: s=%d, D=%d, Xi=%f,  xmin = %f xmax=%f Vmin=%f Vmax=%f\n",s,D,*Xi,  *xmin[0],*xmax[0],*Vmin, *Vmax); 

		
	for (d=0; d<D; d++)
	{
	// 0,1,2  translation
	// 3,4,5,6 quaternion  (not axis angle)
	//  7..   torsion/dihedral in radians

		if(d > 2 && d < 7)
		{
			temp = xmin[d] - xmax[d];
			Vmax[d] = fabs(temp);
			Vmin[d] = -Vmax[d];
            //printf("1: temp = %f: Vmax[%d] = %f\n", fabs(temp), d,Vmax[d]);
		} else if (d > 6)
		{
			Vmax[d] = PI;
			Vmin[d] = -PI;
		} else
		{
			temp = xmin[d] - xmax[d];
			Vmax[d] = fabs(temp) / 2;
			Vmin[d] = -Vmax[d];
		}
		Xi[s].x[d] = random_range(xmin[d], xmax[d]);
        //printf("Xi[%d].x[%d] = %f\n", s,d,Xi[s].x[d]);
		Xi[s].prev_x[d] = Xi[s].x[d];
		Vi[s].v[d] = random_range(Vmin[d], Vmax[d] ); 
        //printf("Vi[%d].v[%d] = %f\n", s,d,Vi[s].v[d]);
	}
    //printf("after loop in initialiseParticle: s=%d, D=%d, Xi=%f,  xmin = %f xmax=%f Vmin=%f Vmax=%f\n",s,D,*Xi,  *xmin[0],*xmax[0],*Vmin, *Vmax); 
    //printf("end initialiseParticle: s=%d, D=%d, Xi=%f  \n",s,D,*Xi);
}
			
void swarmActivity(const int S, const int D, const Position *const Xi, const int nb_eval, const int outlev)
{
	int s, d;
	double swarm_activity;
	double position_dist[PSO_S_MAX];    // Euclidian distance of position and prev position
	double diff;                    // just a helper variable for calculating the distance
	// RG Swarm activity begin
        swarm_activity = 0;
        for(s=0; s<S; s++)
	{
        	for(d=0; d<D; d++)
                	{
				diff = Xi[s].x[d] - Xi[s].prev_x[d];
				//pr(logFile, "swarm_move %d particle= %d dim= %d diff= %f\n", nb_eval, s, d, diff);
				//                        // Differences of angles
				if (d > 5)
				{
					if (diff < -PI) diff = diff + PI;
					if (diff > PI)  diff = diff - PI;
					
				}
				position_dist[s] += (diff * diff);
			}
		position_dist[s] = sqrt(position_dist[s]);
		swarm_activity += position_dist[s];
	}
	if (outlev >= LOGRUNVV )
	{
		swarm_activity = (swarm_activity) / (S * D);
		pr(logFile, "swarm_move %d SA= %f\n", nb_eval+1, swarm_activity);
	}
	// RG Swarm activity end
}

/* EOF */

			
