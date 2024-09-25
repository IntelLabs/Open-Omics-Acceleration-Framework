/*

 $Id: init_gridmap.cpp,v 1.2 2012/04/24 23:33:31 mp Exp $

 AutoGrid 

Copyright (C) 2011 The Scripps Research Institute. All rights reserved.

 AutoGrid is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#include <cstring>
#include <cstdlib>

#include "typedefs.h"
#include "structs.h"
#include "gridmap.h"

MapObject *  init_gridmap( const int num_ligand_types, 
                          const char ligand_types[][3], 
                          const int num_receptor_types, 
                          const char receptor_types[][3], 
                          const int num_maps, 
                          const int num_grid_points_per_map,
                          FILE * logFile
                          )
{
MapObject * gridmap;

        /* num_maps is the number of maps to be created:
         * the number of ligand atom types, plus 1 for the electrostatic map,  if any
         * AutoDock can only read in MAX_MAPS maps, which must include
         * the ligand atom maps and electrostatic map */

        /* Check to see if there is enough memory to store these map objects */
        gridmap = (MapObject *)malloc(sizeof(MapObject) * num_maps);

        if ( gridmap == NULL ) {
            // exits if FATAL_ERROR return code
            print_error( logFile, FATAL_ERROR, "Could not allocate memory to create the MapObject \"gridmap\".\n" );
        }

        // Initialize the gridmap MapObject
        for (int i=0; i<num_maps; i++) {
            gridmap[i].atom_type = 0; /*corresponds to receptor numbers????*/
            gridmap[i].map_index = 0;
            gridmap[i].is_covalent = 0;
            gridmap[i].is_hbonder = 0;
            gridmap[i].map_fileptr = (FILE *)NULL;
            strcpy(gridmap[i].map_filename, "");
            strcpy(gridmap[i].type,""); /*eg HD or OA or NA or N*/
            gridmap[i].constant = 0.0L; /*this will become obsolete*/
            gridmap[i].energy_max = 0.0L;
            gridmap[i].energy_min = 0.0L;
            gridmap[i].energy = 0.0L;
            gridmap[i].vol_probe = 0.0L;
            gridmap[i].solpar_probe = 0.0L;
            gridmap[i].Rij = 0.0L;
            gridmap[i].epsij = 0.0L;
            gridmap[i].hbond = NON; /*hbonding character: */
            gridmap[i].Rij_hb = 0.0L;
            gridmap[i].epsij_hb = 0.0L;
            /*per gridmap[i].receptor type parameters, ordered as in receptor_types*/
            for (int j=0; j<NUM_RECEPTOR_TYPES; j++) {
                gridmap[i].nbp_r[j] = 0.0L; /*radius of energy-well minimum*/
                gridmap[i].nbp_eps[j] = 0.0L;/*depth of energy-well minimum*/
                gridmap[i].xA[j] =0; /*generally 12*/
                gridmap[i].xB[j] =0; /*6 for non-hbonders 10 for h-bonders*/
                gridmap[i].hbonder[j] =0;
            } // j
        } // i


        /* Check to see if the number of grid points requested will be
         * feasible; give warning if not enough memory. */
        if (num_grid_points_per_map > 0) {
            Real *dummy_map;
            dummy_map = (Real *)malloc(sizeof(Real) * (num_maps * num_grid_points_per_map));
            if (!dummy_map) {
                /* Too many maps requested */
                char message[1000];
                (void) sprintf( message, "There may not be enough memory to store these grid maps in AutoDock; \ntry reducing the number of ligand atom types (you have %d including electrostatics) \nor reducing the size of the grid maps; \n or try running AutoDock on a machine with more RAM than this one.\n", num_maps);
                print_error( logFile, WARNING, message );
            } else {
                /* free up this memory right away; we were just testing to
                 * see if we had enough when we try to run AutoDock */
                free(dummy_map);
            }
        } else {
            print_error( logFile, FATAL_ERROR, "You need to set the number of grid points using \"npts\" before setting the ligand atom types, using \"ligand_types\".\n" );
        } /* ZZZZZZZZZZZZZZZZZ*/
        if (!gridmap) {
            char message[1000];
            (void) sprintf( message, "Too many ligand atom types; there is not enough memory to create these maps.  Try using fewer atom types than %d.\n", num_ligand_types);
            print_error( logFile, FATAL_ERROR, message);
        }

        for (int i = 0;  i < num_ligand_types;  i++) {
            ParameterEntry * found_parm;
            gridmap[i].is_covalent = FALSE;
            gridmap[i].is_hbonder = FALSE;
            gridmap[i].map_index = i;
            strcpy(gridmap[i].type, ligand_types[i]); /*eg HD or OA or NA or N*/
            found_parm = apm_find(ligand_types[i]);
            if (strcmp(ligand_types[i],"Z")==0){
                fprintf(logFile, "Found covalent map atomtype\n");
                gridmap[i].is_covalent = TRUE;}
            gridmap[i].atom_type = found_parm->map_index;
            gridmap[i].solpar_probe = found_parm->solpar;
            gridmap[i].vol_probe = found_parm->vol;
            gridmap[i].Rij = found_parm->Rij;
            gridmap[i].epsij = found_parm->epsij;
            gridmap[i].hbond = found_parm->hbond;
            gridmap[i].Rij_hb = found_parm->Rij_hb;
            gridmap[i].epsij_hb = found_parm->epsij_hb;
            if (gridmap[i].hbond>0){       //enum: NON,DS,D1,AS,A1,A2
                gridmap[i].is_hbonder=TRUE;}

#ifdef DEBUG
            (void) fprintf(logFile, " setting ij parms for map %d \n",i);
            (void) fprintf(logFile, "for gridmap[%d], type->%s,Rij->%6.4f, epsij->%6.4f, hbond->%d\n",i,found_parm->autogrid_type, gridmap[i].Rij, gridmap[i].epsij,gridmap[i].hbond);
#endif
            for (int j=0; j<num_receptor_types; j++){
                found_parm = apm_find(receptor_types[j]);
                gridmap[i].nbp_r[j] = (gridmap[i].Rij + found_parm->Rij)/2.;
                gridmap[i].nbp_eps[j] = sqrt(gridmap[i].epsij * found_parm->epsij);
                /*apply the vdW forcefield parameter/weight here */
                // This was removed because "setup_p_l" does this for us... gridmap[i].nbp_eps[j] *= FE_coeff_vdW;
                gridmap[i].xA[j] = 12;
                /*setup hbond dependent stuff*/
                gridmap[i].xB[j] = 6;
                gridmap[i].hbonder[j] = 0;
                if ((int)(gridmap[i].hbond)>2 &&
                        ((int)found_parm->hbond==1||(int)found_parm->hbond==2)){ /*AS,A1,A2 map vs DS,D1 probe*/
                    gridmap[i].xB[j] = 10;
                    gridmap[i].hbonder[j] = 1;
                    //gridmap[i].is_hbonder = TRUE;
                    /*Rij and epsij for this hb interaction in
                     * parm_data.dat file as Rii and epsii for heavy atom
                     * hb factors*/
                    gridmap[i].nbp_r[j] = gridmap[i].Rij_hb;
                    gridmap[i].nbp_eps[j] = gridmap[i].epsij_hb;

                    /*apply the hbond forcefield parameter/weight here */
                    // This was removed because "setup_p_l" does this for us... gridmap[i].nbp_eps[j] *= FE_coeff_hbond;
#ifdef DEBUG
                    (void) fprintf(logFile, "set %d-%d hb eps to %6.4f*%6.4f=%6.4f\n",i,j,gridmap[i].epsij_hb,found_parm->epsij_hb, gridmap[i].nbp_eps[j]);
#endif
                } else if (((int)gridmap[i].hbond==1||(int)gridmap[i].hbond==2) &&
                            ((int)found_parm->hbond>2)) { /*DS,D1 map vs AS,A1,A2 probe*/
                    gridmap[i].xB[j] = 10;
                    gridmap[i].hbonder[j] = 1;
                    gridmap[i].is_hbonder = TRUE;
                    /*Rij and epsij for this hb interaction in
                     * parm_data.dat file as Rii and epsii for heavy atom
                     * hb factors*/
                    gridmap[i].nbp_r[j] = found_parm->Rij_hb;
                    gridmap[i].nbp_eps[j] = found_parm->epsij_hb;

                    /*apply the hbond forcefield parameter/weight here */
                    // This was removed because "setup_p_l" does this for us... gridmap[i].nbp_eps[j] *= FE_coeff_hbond;
#ifdef DEBUG
                    (void) fprintf(logFile, "2: set %d-%d hb eps to %6.4f*%6.4f=%6.4f\n",i,j,gridmap[i].epsij_hb,found_parm->epsij_hb, gridmap[i].nbp_eps[j]);
#endif
                } 
#ifdef DEBUG
                (void) fprintf(logFile, "vs receptor_type[%d]:type->%s, hbond->%d ",j,found_parm->autogrid_type, (int)found_parm->hbond);
                (void) fprintf(logFile, "nbp_r->%6.4f, nbp_eps->%6.4f,xB=%d,hbonder=%d\n",gridmap[i].nbp_r[j], gridmap[i].nbp_eps[j],gridmap[i].xB[j], gridmap[i].hbonder[j]);
#endif
            }; /*initialize energy parms for each possible receptor type*/
        } /*for each map*/
        (void) fprintf( logFile, "\nAtom type names for ligand atom types 1-%d used for ligand-atom affinity grid maps:\n\n", num_ligand_types);
        for (int i = 0;  i < num_ligand_types;  i++) {
            (void) fprintf( logFile, "\t\t\tAtom type number %d corresponds to atom type name \"%s\".\n", gridmap[i].map_index, gridmap[i].type);
            if (gridmap[i].is_covalent == TRUE) {
              (void) fprintf( logFile, "\nAtom type number %d will be used to calculate a covalent affinity grid map\n\n", i + 1);
            }
        }
return gridmap;
}
