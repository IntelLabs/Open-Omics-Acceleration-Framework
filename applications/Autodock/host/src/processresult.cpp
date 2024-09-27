/*

AutoDock-GPU, an OpenCL implementation of AutoDock 4.2 running a Lamarckian Genetic Algorithm
Copyright (C) 2017 TU Darmstadt, Embedded Systems and Applications Group, Germany. All rights reserved.
For some of the code, Copyright (C) 2019 Computational Structural Biology Center, the Scripps Research Institute.

AutoDock is a Trade Mark of the Scripps Research Institute.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/


#include <stdio.h>
#include "processresult.h"


void arrange_result(
                          float* final_population,
                          float* energies,
                    const int    pop_size
                   )
// The function arranges the rows of the input array (first array index is considered to be the row
// index) according to the sum of [] [38] and [][39] elements, which can be used for arranging the
// genotypes of the final population according to the sum of energy values. Genotypes with lower
// energies will be placed at lower row indexes. The second parameter must be equal to the size of
// the population, the arrangement will be performed only on the first pop_size part of final_population.
{
	int i,j;
	float temp_genotype[GENOTYPE_LENGTH_IN_GLOBMEM];
	float temp_energy;

	for (j=0; j<pop_size-1; j++)
		for (i=pop_size-2; i>=j; i--) // arrange according to sum of inter- and intramolecular energies
			if (energies[i] > energies[i+1])
			{
				memcpy(temp_genotype, final_population+i*GENOTYPE_LENGTH_IN_GLOBMEM, GENOTYPE_LENGTH_IN_GLOBMEM*sizeof(float));
				memcpy(final_population+i*GENOTYPE_LENGTH_IN_GLOBMEM, final_population+(i+1)*GENOTYPE_LENGTH_IN_GLOBMEM, GENOTYPE_LENGTH_IN_GLOBMEM*sizeof(float));
				memcpy(final_population+(i+1)*GENOTYPE_LENGTH_IN_GLOBMEM, temp_genotype, GENOTYPE_LENGTH_IN_GLOBMEM*sizeof(float));

				temp_energy = energies[i];
				energies[i] = energies[i+1];
				energies[i+1] = temp_energy;
			}
}


void write_basic_info(
                            FILE*       fp,
                      const Liganddata* ligand_ref,
                      const Dockpars*   mypars,
                      const Gridinfo*   mygrid,
                      const int*        argc,
                            char**      argv
                     )
// The function writes basic information (such as docking parameters) to the file whose file pointer is the first parameter of the function.
{
	int i;

	fprintf(fp, "***********************************\n");
	fprintf(fp, "**    AUTODOCK-GPU REPORT FILE   **\n");
	fprintf(fp, "***********************************\n\n\n");

	// Writing out docking parameters

	fprintf(fp, "         DOCKING PARAMETERS        \n");
	fprintf(fp, "===================================\n\n");

	fprintf(fp, "Ligand file:                               %s\n", mypars->ligandfile);
	bool flexres = false;
	if (mypars->flexresfile!=NULL){
			if ( strlen(mypars->flexresfile)>0 ) {
				fprintf(fp, "Flexible residue file:                     %s\n", mypars->flexresfile);
				flexres = true;
			}
	}
	fprintf(fp, "Grid fld file:                             %s\n", mypars->fldfile);

	fprintf(fp, "Number of energy evaluations:              %ld\n", mypars->num_of_energy_evals);
	fprintf(fp, "Number of generations:                     %ld\n", mypars->num_of_generations);
	fprintf(fp, "Size of population:                        %ld\n", mypars->pop_size);
	fprintf(fp, "Rate of crossover:                         %lf%%\n", (double) mypars->crossover_rate);
	fprintf(fp, "Tournament selection probability limit:    %lf%%\n", (double) mypars->tournament_rate);
	fprintf(fp, "Rate of mutation:                          %lf%%\n", (double) mypars->mutation_rate);
	fprintf(fp, "Maximal allowed delta movement:            +/- %lfA\n", (double) mypars->abs_max_dmov*mygrid->spacing);
	fprintf(fp, "Maximal allowed delta angle:               +/- %lf\n\n", (double) mypars->abs_max_dang);

	fprintf(fp, "Rate of local search:                      %lf%%\n", mypars->lsearch_rate);

	fprintf(fp, "Maximal number of local search iterations: %ld\n", mypars->max_num_of_iters);
	fprintf(fp, "Rho lower bound:                           %lf\n", (double) mypars->rho_lower_bound);
	fprintf(fp, "Spread of local search delta movement:     %lfA\n", (double) mypars->base_dmov_mul_sqrt3*mygrid->spacing/sqrt(3.0));
	fprintf(fp, "Spread of local search delta angle:        %lf\n", (double) mypars->base_dang_mul_sqrt3/sqrt(3.0));
	fprintf(fp, "Limit of consecutive successes/failures:   %ld\n\n", mypars->cons_limit);

	fprintf(fp, "Unbound model:                             ");
	if (mypars->unbound_model == 0)
		fprintf(fp, "BOUND\n");
	else
		if (mypars->unbound_model == 1)
			fprintf(fp, "EXTENDED\n");
		else
			fprintf(fp, "COMPACT\n");

	fprintf(fp, "Number of pdb files to be generated:       %d\n", mypars->gen_pdbs);

	fprintf(fp, "Initial population:                        ");
	if (!mypars->load_xml)
		fprintf(fp, "GENERATE\n");
	else
		fprintf(fp, "LOAD FROM FILE (%s)\n",mypars->load_xml);

	fprintf(fp, "\n\nProgram call in command line was:          ");
	for (i=0; i<*argc; i++)
		fprintf(fp, "%s ", argv [i]);
	fprintf(fp, "\n\n\n");

	// Writing out receptor parameters

	fprintf(fp, "        RECEPTOR PARAMETERS        \n");
	fprintf(fp, "===================================\n\n");

	fprintf(fp, "Receptor name:                             %s\n", mygrid->receptor_name);
	fprintf(fp, "Number of grid points (x, y, z):           %d, %d, %d\n", mygrid->size_xyz [0], mygrid->size_xyz [1], mygrid->size_xyz [2]);
	fprintf(fp, "Grid size (x, y, z):                       %lf, %lf, %lfA\n", mygrid->size_xyz_angstr [0], mygrid->size_xyz_angstr [1], mygrid->size_xyz_angstr [2]);
	fprintf(fp, "Grid spacing:                              %lfA\n", mygrid->spacing);
	fprintf(fp, "\n\n");

	// Writing out ligand parameters
	if(flexres)
		fprintf(fp, "     LIGAND+FLEXRES PARAMETERS     \n");
	else
		fprintf(fp, "         LIGAND PARAMETERS         \n");
	fprintf(fp, "===================================\n\n");

	fprintf(fp, "Ligand name:                               ");
	int len = strlen(mypars->ligandfile) - 6;
	for(i=0; i<len; i++) fputc(mypars->ligandfile[i], fp);
	fputc('\n', fp);
	if(flexres){
		fprintf(fp, "Flexres name:                              ");
		int len = strlen(mypars->flexresfile) - 6;
		for(i=0; i<len; i++) fputc(mypars->flexresfile[i], fp);
		fputc('\n', fp);
		fprintf(fp, "Number of ligand atoms:                    %d\n", ligand_ref->true_ligand_atoms);
		fprintf(fp, "Number of flexres atoms:                   %d\n", ligand_ref->num_of_atoms-ligand_ref->true_ligand_atoms);
		fprintf(fp, "Number of ligand rotatable bonds:          %d\n", ligand_ref->true_ligand_rotbonds);
		fprintf(fp, "Number of flexres rotatable bonds:         %d\n", ligand_ref->num_of_rotbonds-ligand_ref->true_ligand_rotbonds);
	}
	fprintf(fp, "Number of atoms:                           %d\n", ligand_ref->num_of_atoms);
	fprintf(fp, "Number of rotatable bonds:                 %d\n", ligand_ref->num_of_rotbonds);
	fprintf(fp, "Number of atom types:                      %d\n", ligand_ref->num_of_atypes);

	fprintf(fp, "Number of intraE contributors:             %d\n", ligand_ref->num_of_intraE_contributors);
	fprintf(fp, "Number of required rotations:              %d\n", ligand_ref->num_of_rotations_required);
	fprintf(fp, "Number of rotation cycles:                 %d\n", ligand_ref->num_of_rotcyc);

	fprintf(fp, "\n\n");
}

void write_basic_info_dlg(
                                FILE*       fp,
                          const Liganddata* ligand_ref,
                          const Dockpars*   mypars,
                          const Gridinfo*   mygrid,
                          const int*        argc,
                                char**      argv
                         )
// The function writes basic information (such as docking parameters) to the file whose file pointer is the first parameter of the function.
{
	int i;

	if(mypars->xml2dlg && mypars->dlg2stdout) fprintf(fp, "\nXML2DLG: %s\n", mypars->load_xml);
	fprintf(fp, "AutoDock-GPU version: %s\n\n", VERSION);

	fprintf(fp, "**********************************************************\n");
	fprintf(fp, "**    AutoDock-GPU AUTODOCKTOOLS-COMPATIBLE DLG FILE    **\n");
	fprintf(fp, "**********************************************************\n\n\n");

	// Writing out docking parameters

	fprintf(fp, "    DOCKING PARAMETERS\n");
	fprintf(fp, "    ________________________\n\n\n");

	fprintf(fp, "Ligand file:                               %s\n", mypars->ligandfile);
	bool flexres = false;
	if (mypars->flexresfile!=NULL){
			if ( strlen(mypars->flexresfile)>0 ) {
				fprintf(fp, "Flexible residue file:                     %s\n", mypars->flexresfile);
				flexres = true;
			}
	}
	fprintf(fp, "Grid fld file:                             %s\n\n", mypars->fldfile);

	fprintf(fp, "Random seed:                               %u", mypars->seed[0]);
	if(mypars->seed[1]>0) fprintf(fp,", %u",mypars->seed[1]);
	if(mypars->seed[2]>0) fprintf(fp,", %u",mypars->seed[2]);
	fprintf(fp, "\n");
	fprintf(fp, "Number of runs:                            %lu\n", mypars->num_of_runs);
	
	if(!mypars->xml2dlg){
		fprintf(fp, "Number of energy evaluations:              %ld\n", mypars->num_of_energy_evals);
		fprintf(fp, "Number of generations:                     %ld\n", mypars->num_of_generations);
		fprintf(fp, "Size of population:                        %ld\n", mypars->pop_size);
		fprintf(fp, "Rate of crossover:                         %lf%%\n", (double) mypars->crossover_rate);
		fprintf(fp, "Tournament selection probability limit:    %lf%%\n", (double) mypars->tournament_rate);
		fprintf(fp, "Rate of mutation:                          %lf%%\n", (double) mypars->mutation_rate);
		fprintf(fp, "Maximal allowed delta movement:            +/- %lfA\n", (double) mypars->abs_max_dmov*mygrid->spacing);
		fprintf(fp, "Maximal allowed delta angle:               +/- %lf\n\n", (double) mypars->abs_max_dang);

		fprintf(fp, "Rate of local search:                      %lf%%\n", mypars->lsearch_rate);

		fprintf(fp, "Maximal number of local search iterations: %ld\n", mypars->max_num_of_iters);
		fprintf(fp, "Rho lower bound:                           %lf\n", (double) mypars->rho_lower_bound);
		fprintf(fp, "Spread of local search delta movement:     %lfA\n", (double) mypars->base_dmov_mul_sqrt3*mygrid->spacing/sqrt(3.0));
		fprintf(fp, "Spread of local search delta angle:        %lf\n", (double) mypars->base_dang_mul_sqrt3/sqrt(3.0));
		fprintf(fp, "Limit of consecutive successes/failures:   %ld\n\n", mypars->cons_limit);
	}

		fprintf(fp, "Handle symmetry during clustering:         ");
	if (mypars->handle_symmetry)
		fprintf(fp, "YES\n");
	else
		fprintf(fp, "NO\n");

	fprintf(fp, "RMSD tolerance:                            %lfA\n\n", mypars->rmsd_tolerance);

	if(!mypars->xml2dlg){ // This is necessary to avoid excruciatingly long command line outputs with wild cards (like *.xml)
		fprintf(fp, "Program call in command line was:          ");
		for (i=0; i<*argc; i++)
			fprintf(fp, "%s ", argv [i]);
		fprintf(fp, "\n\n\n");
	}

	// Writing out receptor parameters

	fprintf(fp, "    GRID PARAMETERS\n");
	fprintf(fp, "    ________________________\n\n\n");

	fprintf(fp, "Receptor name:                             %s\n", mygrid->receptor_name);
	fprintf(fp, "Number of grid points (x, y, z):           %d, %d, %d\n", mygrid->size_xyz [0],
			mygrid->size_xyz [1], mygrid->size_xyz [2]);
	fprintf(fp, "Grid size (x, y, z):                       %lf, %lf, %lfA\n", mygrid->size_xyz_angstr [0],
			mygrid->size_xyz_angstr [1], mygrid->size_xyz_angstr [2]);
	fprintf(fp, "Grid spacing:                              %lfA\n", mygrid->spacing);
	fprintf(fp, "\n\n");

	// Writing out ligand parameters
	if(flexres)
	fprintf(fp, "    LIGAND+FLEXRES PARAMETERS\n");
	else
		fprintf(fp, "    LIGAND PARAMETERS\n");
	fprintf(fp, "    ________________________\n\n\n");

	fprintf(fp, "Ligand name:                               ");
	int len = strlen(mypars->ligandfile) - 6;
	for(i=0; i<len; i++) fputc(mypars->ligandfile[i], fp);
	fputc('\n', fp);
	if(flexres){
		fprintf(fp, "Flexres name:                              ");
		int len = strlen(mypars->flexresfile) - 6;
		for(i=0; i<len; i++) fputc(mypars->flexresfile[i], fp);
		fputc('\n', fp);
		fprintf(fp, "Number of ligand atoms:                    %d\n", ligand_ref->true_ligand_atoms);
		fprintf(fp, "Number of flexres atoms:                   %d\n", ligand_ref->num_of_atoms-ligand_ref->true_ligand_atoms);
		fprintf(fp, "Number of ligand rotatable bonds:          %d\n", ligand_ref->true_ligand_rotbonds);
		fprintf(fp, "Number of flexres rotatable bonds:         %d\n", ligand_ref->num_of_rotbonds-ligand_ref->true_ligand_rotbonds);
	}
	fprintf(fp, "Number of atoms:                           %d\n", ligand_ref->num_of_atoms);
	fprintf(fp, "Number of rotatable bonds:                 %d\n", ligand_ref->num_of_rotbonds);
	fprintf(fp, "Number of atom types:                      %d\n", ligand_ref->num_of_atypes);
	fprintf(fp, "\n\n");

	if(!mypars->xml2dlg){
		fprintf(fp, "    DUMMY DATA (only for ADT-compatibility)\n");
		fprintf(fp, "    ________________________\n\n\n");
		fprintf(fp, "DPF> outlev 1\n");
		fprintf(fp, "DPF> ga_run %lu\n", mypars->num_of_runs);
		fprintf(fp, "DPF> fld %s.maps.fld\n", mygrid->receptor_name);
		fprintf(fp, "DPF> move %s\n", mypars->ligandfile);
		if(flexres) fprintf(fp, "DPF> flexres %s\n", mypars->flexresfile);
		fprintf(fp, "\n\n");
	}
}

void make_resfiles(
                         float*        final_population,
                         float*        energies,
                   const Liganddata*   ligand_ref,
                   const Liganddata*   ligand_from_pdb,
                   const Liganddata*   ligand_xray,
                   const Dockpars*     mypars,
                         int           evals_performed,
                         int           generations_used,
                   const Gridinfo*     mygrid,
                   const float*        grids,
                   const int*          argc,
                         char**        argv,
                         int           debug,
                         int           run_cnt,
                         float&        best_energy_of_all,
                         Ligandresult* best_result
                  )
// The function writes out final_population generated by get_result
// as well as different parameters about the docking, the receptor and the ligand to a file called fdock_report.txt in a
// readable and understandable format. The ligand_from_pdb parametere must be the Liganddata which includes the original
// ligand conformation as the result conformations will be compared to this one. The structs containing the grid informations
// and docking parameters are requrided as well as the number and values of command line arguments. The ligand_ref parameter
// describes the ligand with the reference orientation (gene values of final_population refer to this one, that is, this can
// be moved and rotated according to the genotype values). The function returns some information about the best result wich
// was found with the best_result parameter.
{
	FILE* fp = stdout; // takes care of compile warning down below (and serves as a visual bug tracker in case fp is written to accidentally)
	int i,j;
	double entity_rmsds;
	Liganddata temp_docked;
	int len = strlen(mypars->ligandfile) - 6 + 24 + 10 + 10; // length with added bits for things below (numbers below 11 digits should be a safe enough threshold)
	char* temp_filename = (char*)malloc((len+1)*sizeof(char)); // +\0 at the end
	char* name_ext_start;
	std::vector<AnalysisData> analysis;
	float accurate_interE;
	float accurate_intraflexE;
	float accurate_intraE;
	float accurate_interflexE;

	int pop_size = mypars->pop_size;

	sprintf(temp_filename, "final_population_run%d.txt", run_cnt+1);

	if (mypars->gen_finalpop) // if final population files are not required, no file will be opened.
	{
		fp = fopen(temp_filename, "w");

		write_basic_info(fp, ligand_ref, mypars, mygrid, argc, argv); // Write basic information about docking and molecule parameters to file

		fprintf(fp, "           COUNTER STATES           \n");
		fprintf(fp, "===================================\n\n");
		fprintf(fp, "Number of energy evaluations performed:    %d\n", evals_performed);
		fprintf(fp, "Number of generations used:                %d\n", generations_used);
		fprintf(fp, "\n\n");
		fprintf(fp, "     STATE OF FINAL POPULATION     \n");
		fprintf(fp, "===================================\n\n");

		fprintf(fp, " Entity |      dx [A]      |      dy [A]      |      dz [A]      |      phi []      |     theta []     |  alpha_genrot [] |");
		for (i=0; i<ligand_from_pdb->num_of_rotbonds; i++)
			fprintf(fp, "  alpha_rotb%2d [] |", i);
		fprintf(fp, " intramolecular energy | intermolecular energy |     total energy calculated by CPU / calculated by GPU / difference    | RMSD [A] | \n");

		fprintf(fp, "--------+------------------+------------------+------------------+------------------+------------------+------------------+");
		for (i=0; i<ligand_from_pdb->num_of_rotbonds; i++)
			fprintf(fp, "------------------+");
		fprintf(fp, "-----------------------+-----------------------+------------------------------------------------------------------------+----------+ \n");
	}

	// Writing out state of final population

	strcpy(temp_filename, mypars->ligandfile);
	name_ext_start = temp_filename + strlen(mypars->ligandfile) - 6; // without .pdbqt

	IntraTables tables(ligand_ref, mypars->coeffs.scaled_AD4_coeff_elec, mypars->coeffs.AD4_coeff_desolv, mypars->qasp);
	for (i=0; i<pop_size; i++)
	{
		temp_docked = *ligand_ref;
		
		if(mypars->xml2dlg){
			double axisangle[4];
			double genotype [GENOTYPE_LENGTH_IN_GLOBMEM];
			for (unsigned int j=0; j<ACTUAL_GENOTYPE_LENGTH; j++)
				genotype [j] = (final_population+i*GENOTYPE_LENGTH_IN_GLOBMEM)[j];
			genotype[GENOTYPE_LENGTH_IN_GLOBMEM-1] = (final_population+i*GENOTYPE_LENGTH_IN_GLOBMEM)[GENOTYPE_LENGTH_IN_GLOBMEM-1];
			axisangle[0] = genotype[3];
			axisangle[1] = genotype[4];
			axisangle[2] = genotype[5];
			axisangle[3] = genotype[GENOTYPE_LENGTH_IN_GLOBMEM-1];
			change_conform(&temp_docked, mygrid, genotype, axisangle, debug);
		} else{
			change_conform_f(&temp_docked, mygrid, final_population+i*GENOTYPE_LENGTH_IN_GLOBMEM, debug);
		}
		
		// the map interaction of flex res atoms is stored in accurate_intraflexE
		accurate_interE = calc_interE_f(mygrid, &temp_docked, grids, 0.0005, debug, accurate_intraflexE);	//calculating the intermolecular energy

		if (mypars->contact_analysis && (i==0)){
			analysis = analyze_ligand_receptor(mygrid, &temp_docked, mypars->receptor_atoms.data(), mypars->receptor_map, mypars->receptor_map_list, 0.0005, debug, mypars->H_cutoff, mypars->V_cutoff);
		}

		if (i == 0) // additional calculations for ADT-compatible result file, only in case of best conformation
			calc_interE_peratom_f(mygrid, &temp_docked, grids, 0.0005, &(best_result->interE_elec), best_result->peratom_vdw, best_result->peratom_elec, debug);

		scale_ligand(&temp_docked, mygrid->spacing);
		
		// the interaction between flex res and ligand is stored in accurate_interflexE
		if(mypars->contact_analysis && (i==0))
			accurate_intraE = calc_intraE_f(&temp_docked, 8, mypars->smooth, 0, mypars->elec_min_distance, tables, debug, accurate_interflexE, mypars->nr_mod_atype_pairs, mypars->mod_atype_pairs, &analysis, mypars->receptor_atoms.data() + mypars->nr_receptor_atoms, mypars->R_cutoff, mypars->H_cutoff, mypars->V_cutoff);
		else
			accurate_intraE = calc_intraE_f(&temp_docked, 8, mypars->smooth, 0, mypars->elec_min_distance, tables, debug, accurate_interflexE, mypars->nr_mod_atype_pairs, mypars->mod_atype_pairs);

		move_ligand(&temp_docked, mygrid->origo_real_xyz, mygrid->origo_real_xyz); //moving it according to grid location

//		for (unsigned int atom_id=0; atom_id < temp_docked.num_of_atoms; atom_id++)
//			printf("%i: %lf, %lf, %lf\n", atom_id+1, temp_docked.atom_idxyzq [atom_id][1], temp_docked.atom_idxyzq [atom_id][2], temp_docked.atom_idxyzq [atom_id][3]);

		if (mypars->given_xrayligandfile == true) {
			entity_rmsds = calc_rmsd(ligand_xray, &temp_docked, mypars->handle_symmetry); //calculating rmds compared to original xray file
		}
		else {
			entity_rmsds = calc_rmsd(ligand_from_pdb, &temp_docked, mypars->handle_symmetry); //calculating rmds compared to original pdb file
		}

		// copying best result to output parameter
		if (i == 0) // assuming this is the best one (final_population is arranged), however, the
		{           // arrangement was made according to the inaccurate values calculated by FPGA
			best_result->genotype = final_population+i*GENOTYPE_LENGTH_IN_GLOBMEM;
			best_result->interE = accurate_interE;
			best_result->interflexE = accurate_interflexE;
			best_result->intraE = accurate_intraE;
			best_result->intraflexE = accurate_intraflexE;
			best_result->reslig_realcoord = temp_docked;
			best_result->rmsd_from_ref = entity_rmsds;
			best_result->run_number = run_cnt+1;
			if(mypars->contact_analysis) best_result->analysis = analysis;
		}

		// generating best.pdbqt
		if (i == 0)
			if (best_energy_of_all > accurate_interE + accurate_intraE)
			{
				best_energy_of_all = accurate_interE + accurate_intraE;

				if (mypars->gen_best)
					gen_new_pdbfile(mypars->ligandfile, "best.pdbqt", &temp_docked);
			}

		if (i < mypars->gen_pdbs) //if it is necessary, making new pdbqts for best entities
		{
			sprintf(name_ext_start, "_docked_run%d_entity%d.pdbqt", run_cnt+1, i+1); //name will be <original pdb filename>_docked_<number starting from 1>.pdb
			gen_new_pdbfile(mypars->ligandfile, temp_filename, &temp_docked);
		}
		if (mypars->gen_finalpop)
		{
			fprintf(fp, "  %3d   |", i+1);

			for (j=0; j<3; j++)
				fprintf(fp, "    %10.3f    |", final_population [i*GENOTYPE_LENGTH_IN_GLOBMEM+j]*(mygrid->spacing));
			for (j=3; j<6+ligand_from_pdb->num_of_rotbonds; j++)
				fprintf(fp, "    %10.3f    |", final_population [i*GENOTYPE_LENGTH_IN_GLOBMEM+j]);

			fprintf(fp, " %21.3f |", accurate_intraE);
			fprintf(fp, " %21.3f |", accurate_interE);
			fprintf(fp, "  %21.3f / %21.3f / %21.3f |", accurate_intraE + accurate_interE, energies[i], energies[i] - (accurate_intraE + accurate_interE));

			fprintf(fp, " %8.3lf | \n", entity_rmsds);
		}
	}
	if (mypars->gen_finalpop) fclose(fp);
	free(temp_filename);
}

void cluster_analysis(
                            Ligandresult myresults [],
                            int          num_of_runs,
                            char*        report_file_name,
                      const Liganddata*  ligand_ref,
                      const Dockpars*    mypars,
                      const Gridinfo*    mygrid,
                      const int*         argc,
                            char**       argv,
                      const double       docking_avg_runtime,
                      const double program_runtime
                     )
// The function performs ranked cluster analisys similar to that of AutoDock and creates a file with report_file_name name, the result
// will be written to it.
{
	int i,j;
	Ligandresult temp_ligres;
	int num_of_clusters;
	int current_clust_center;
	double temp_rmsd;
	double cluster_tolerance = 2;
	int result_clustered;
	int subrank;
	FILE* fp;
	int cluster_sizes [1000];
	double sum_energy [1000];
	double best_energy [1000];

	const double AD4_coeff_tors = mypars->coeffs.AD4_coeff_tors;
	double torsional_energy;

	// first of all, let's calculate the constant torsional free energy term
	torsional_energy = AD4_coeff_tors * ligand_ref->true_ligand_rotbonds;

	// arranging results according to energy, myresults [0] will be the best one (with lowest energy)
	for (j=0; j<num_of_runs-1; j++)
		for (i=num_of_runs-2; i>=j; i--) // arrange according to sum of inter- and intramolecular energies
			if ((myresults [i]).interE /*+ (myresults [i]).intraE*/ > (myresults [i+1]).interE /*+ (myresults [i+1]).intraE*/) // mimics the behaviour of AD4 unbound_same_as_bound
			//if ((myresults [i]).interE + (myresults [i]).intraE > (myresults [i+1]).interE + (myresults [i+1]).intraE)
			{
				temp_ligres = myresults [i];
				myresults [i] = myresults [i+1];
				myresults [i+1] = temp_ligres;
			}

	for (i=0; i<num_of_runs; i++)
	{
		(myresults [i]).clus_id = 0; // indicates that it hasn't been put into cluster yet
	}

	// the best result is the center of the first cluster
	(myresults [0]).clus_id = 1;
	(myresults [0]).rmsd_from_cluscent = 0;
	num_of_clusters = 1;

	for (i=1; i<num_of_runs; i++) // for each result
	{
		current_clust_center = 0;
		result_clustered = 0;
		for (j=0; j<i; j++) // results with lower id-s are clustered, look for cluster centers
		{
			if ((myresults [j]).clus_id > current_clust_center) // it is the center of a new cluster
			{
				current_clust_center = (myresults [j]).clus_id;
				temp_rmsd = calc_rmsd(&((myresults [j]).reslig_realcoord), &((myresults [i]).reslig_realcoord), mypars->handle_symmetry); // comparing current result with cluster center
				if (temp_rmsd <= cluster_tolerance) // in this case we put result i to cluster with center j
				{
					(myresults [i]).clus_id = current_clust_center;
					(myresults [i]).rmsd_from_cluscent = temp_rmsd;
					result_clustered = 1;
					break;
				}
			}
		}
		if (result_clustered != 1) // if no suitable cluster was found, this is the center of a new one
		{
			num_of_clusters++;
			(myresults [i]).clus_id = num_of_clusters; // new cluster id
			(myresults [i]).rmsd_from_cluscent = 0;
		}
	}

	for (i=1; i<=num_of_clusters; i++) // printing cluster info to file
	{
		subrank = 0;
		cluster_sizes [i-1] = 0;
		sum_energy [i-1] = 0;
		for (j=0; j<num_of_runs; j++)
			if (myresults [j].clus_id == i)
			{
				subrank++;
				(cluster_sizes [i-1])++;
				sum_energy [i-1] += (myresults [j]).interE + /*(myresults [j]).intraE +*/ torsional_energy; // intraE can be commented when unbound_same_as_bound
				(myresults [j]).clus_subrank = subrank;
				if (subrank == 1)
					best_energy [i-1] = (myresults [j]).interE + /*(myresults [j]).intraE +*/ torsional_energy; // intraE can be commented when unbound_same_as_bound
			}
	}

	fp = fopen(report_file_name, "w");

	write_basic_info(fp, ligand_ref, mypars, mygrid, argc, argv); // Write basic information about docking and molecule parameters to file
	fprintf(fp, "           RUN TIME INFO           \n");
	fprintf(fp, "===================================\n\n");

	fprintf(fp, "Average GPU run time for 1 run:           %lfs\n", docking_avg_runtime);
	fprintf(fp, "Total GPU docking run time:               %fs\n", docking_avg_runtime*mypars->num_of_runs);

	fprintf(fp, "Program run time:                          %lfs\n", program_runtime);
	fprintf(fp, "\n\n");

	fprintf(fp, "       CLUSTERING HISTOGRAM        \n");
	fprintf(fp, "===================================\n\n");
	fprintf(fp, " Cluster rank | Num in cluster |   Best energy   |   Mean energy   |    5    10   15   20   25   30   35\n");
	fprintf(fp, "--------------+----------------+-----------------+-----------------+----+----+----+----+----+----+----+\n");

	for (i=1; i<=num_of_clusters; i++)
	{
		fprintf(fp, "      %3d     |       %3d      | %15.3lf | %15.3lf |", i, cluster_sizes [i-1], best_energy [i-1], sum_energy [i-1]/cluster_sizes [i-1]);

		for (j=0; j<cluster_sizes [i-1]; j++)
			fprintf(fp, "#");

		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n");

	fprintf(fp, "              CLUSTERS             \n");
	fprintf(fp, "===================================\n\n");
	fprintf(fp, " Rank | Subrank | Run | Intermolecular E | Intramolecular E | Torsional energy |   Total energy   | Cluster RMSD | Reference RMSD |\n");
	fprintf(fp, "------+---------+-----+------------------+------------------+------------------+------------------+--------------+----------------+\n");

	for (i=1; i<=num_of_clusters; i++) // printing cluster info to file
	{
		for (j=0; j<num_of_runs; j++)
			if (myresults [j].clus_id == i)
			{
				fprintf(fp, "  %3d |   %3d   | %3d |  %15.3lf |  %15.3lf |  %15.3lf |  %15.3lf |     %4.2lf     |      %4.2lf      |\n", (myresults [j]).clus_id, (myresults [j]).clus_subrank, (myresults [j]).run_number,
				            (myresults [j]).interE, (myresults [j]).intraE, torsional_energy, (myresults [j]).interE + /*(myresults [j]).intraE +*/ torsional_energy, (myresults [j]).rmsd_from_cluscent, (myresults [j]).rmsd_from_ref); // intraE can be commented when unbound_same_as_bound
			}
	}
	fclose(fp);
}

void clusanal_gendlg(
                           Ligandresult  myresults [],
                           int           num_of_runs,
                     const Liganddata*   ligand_ref,
                     const Dockpars*     mypars,
                     const Gridinfo*     mygrid,
                     const int*          argc,
                           char**        argv,
                     const double        docking_avg_runtime,
                           unsigned long generations_used,
                           unsigned long evals_performed,
                           double        exec_time,
                           double        idle_time
                    )
// The function performs ranked cluster analysis similar to that of AutoDock and creates a file with report_file_name name, the result
// will be written to it.
{
	int i, j, atom_cnt;
	Ligandresult temp_ligres;
	int num_of_clusters;
	int current_clust_center;
	double temp_rmsd;
	int result_clustered;
	int subrank;
	FILE* fp = stdout;
	FILE* fp_orig;
	FILE* fp_xml;
	int cluster_sizes [1000];
	double sum_energy [1000];
	double best_energy [1000];
	int best_energy_runid [1000];
	char tempstr [256];

	double cluster_tolerance = mypars->rmsd_tolerance;
	const double AD4_coeff_tors = mypars->coeffs.AD4_coeff_tors;
	double torsional_energy;

	// first of all, let's calculate the constant torsional free energy term
	torsional_energy = AD4_coeff_tors * ligand_ref->true_ligand_rotbonds;

	int len = strlen(mypars->resname) + 4 + 1;
	// GENERATING DLG FILE
	if(mypars->output_dlg){
		if(!mypars->dlg2stdout){
			char* report_file_name = (char*)malloc(len*sizeof(char));
			strcpy(report_file_name, mypars->resname);
			strcat(report_file_name, ".dlg");
			fp = fopen(report_file_name, "w");
			free(report_file_name);
		}

		// writing basic info
		write_basic_info_dlg(fp, ligand_ref, mypars, mygrid, argc, argv);

		if(!mypars->xml2dlg){
			fprintf(fp, "           COUNTER STATES           \n");
			fprintf(fp, "___________________________________\n\n");
			fprintf(fp, "Number of energy evaluations performed:    %lu\n", evals_performed);
			fprintf(fp, "Number of generations used:                %lu\n", generations_used);
			fprintf(fp, "\n\n");
		}

		std::string pdbqt_template;
		std::vector<unsigned int> atom_data;
		char lineout [264];
		// writing input pdbqt file
		fprintf(fp, "    INPUT LIGAND PDBQT FILE:\n    ________________________\n\n\n");
		fp_orig = fopen(mypars->ligandfile, "rb"); // fp_orig = fopen(mypars->ligandfile, "r");
		while (fgets(tempstr, 255, fp_orig) != NULL) // reading original ligand pdb line by line
		{
			fprintf(fp, "INPUT-LIGAND-PDBQT: %s", tempstr);
			if ((strncmp("ATOM", tempstr, 4) == 0) || (strncmp("HETATM", tempstr, 6) == 0))
			{
				tempstr[30] = '\0';
				sprintf(lineout, "DOCKED: %s", tempstr);
				pdbqt_template += lineout;
				atom_data.push_back(pdbqt_template.length());
			} else{
				if (strncmp("ROOT", tempstr, 4) == 0)
				{
					pdbqt_template += "DOCKED: USER                              x       y       z     vdW  Elec       q    Type\n";
					pdbqt_template += "DOCKED: USER                           _______ _______ _______ _____ _____    ______ ____\n";
				}
				sprintf(lineout, "DOCKED: %s", tempstr);
				pdbqt_template += lineout;
			}
		}
		fprintf(fp, "\n\n");
		fclose(fp_orig);
		
		// writing input flexres pdbqt file if specified
		if (mypars->flexresfile!=NULL) {
			if ( strlen(mypars->flexresfile)>0 ) {
				fprintf(fp, "    INPUT FLEXRES PDBQT FILE:\n    ________________________\n\n\n");
				fp_orig = fopen(mypars->flexresfile, "rb"); // fp_orig = fopen(mypars->flexresfile, "r");
				while (fgets(tempstr, 255, fp_orig) != NULL) // reading original flexres pdb line by line
				{
					fprintf(fp, "INPUT-FLEXRES-PDBQT: %s", tempstr);
					if ((strncmp("ATOM", tempstr, 4) == 0) || (strncmp("HETATM", tempstr, 6) == 0))
					{
						tempstr[30] = '\0';
						sprintf(lineout, "DOCKED: %s", tempstr);
						pdbqt_template += lineout;
						atom_data.push_back(pdbqt_template.length());
					} else{
						if (strncmp("ROOT", tempstr, 4) == 0)
						{
							pdbqt_template += "DOCKED: USER                              x       y       z     vdW  Elec       q    Type\n";
							pdbqt_template += "DOCKED: USER                           _______ _______ _______ _____ _____    ______ ____\n";
							}
						sprintf(lineout, "DOCKED: %s", tempstr);
						pdbqt_template += lineout;
					}
				}
				fprintf(fp, "\n\n");
				fclose(fp_orig);
			}
		}
		
		// writing docked conformations
		std::string curr_model;
		for (i=0; i<num_of_runs; i++)
		{
			fprintf(fp, "    FINAL DOCKED STATE:\n    ________________________\n\n\n");

			fprintf(fp, "Run:   %d / %lu\n", i+1, mypars->num_of_runs);
			fprintf(fp, "Time taken for this run:   %.3lfs\n\n", docking_avg_runtime);

			if(mypars->contact_analysis){
				if(myresults[i].analysis.size()>0){
					// sort by analysis type
					AnalysisData temp;
					for(unsigned int j=0; j<myresults[i].analysis.size(); j++)
						for(unsigned int k=0; k<myresults[i].analysis.size()-j-1; k++)
							if(myresults[i].analysis[k].type>myresults[i].analysis[k+1].type){
								temp = myresults[i].analysis[k];
								myresults[i].analysis[k]   = myresults[i].analysis[k+1];
								myresults[i].analysis[k+1] = temp;
							}
					fprintf(fp, "ANALYSIS: COUNT %lu\n", myresults[i].analysis.size());
					std::string types    = "TYPE    {";
					std::string lig_id   = "LIGID   {";
					std::string ligname  = "LIGNAME {";
					std::string rec_id   = "RECID   {";
					std::string rec_name = "RECNAME {";
					std::string residue  = "RESIDUE {";
					std::string res_id   = "RESID   {";
					std::string chain    = "CHAIN   {";
					char item[8], pad[8];
					for(unsigned int j=0; j<myresults[i].analysis.size(); j++){
						if(j>0){
							types    += ",";
							lig_id   += ",";
							ligname  += ",";
							rec_id   += ",";
							rec_name += ",";
							residue  += ",";
							res_id   += ",";
							chain    += ",";
						}
						switch(myresults[i].analysis[j].type){
							case 0: types += "   \"R\"";
							        break;
							case 1: types += "   \"H\"";
							        break;
							default:
							case 2: types += "   \"V\"";
							        break;
						}
						sprintf(item, "%5d ", myresults[i].analysis[j].lig_id);   lig_id+=item;
						sprintf(item, "\"%s\"", myresults[i].analysis[j].lig_name); sprintf(pad, "%6s", item); ligname+=pad;
						sprintf(item, "%5d ", myresults[i].analysis[j].rec_id);   rec_id+=item;
						sprintf(item, "\"%s\"", myresults[i].analysis[j].rec_name); sprintf(pad, "%6s", item); rec_name+=pad;
						sprintf(item, "\"%s\"", myresults[i].analysis[j].residue); sprintf(pad, "%6s", item);  residue+=pad;
						sprintf(item, "%5d ", myresults[i].analysis[j].res_id);   res_id+=item;
						sprintf(item, "\"%s\"", myresults[i].analysis[j].chain); sprintf(pad, "%6s", item);    chain+=pad;
					}
					fprintf(fp, "ANALYSIS: %s}\n", types.c_str());
					fprintf(fp, "ANALYSIS: %s}\n", lig_id.c_str());
					fprintf(fp, "ANALYSIS: %s}\n", ligname.c_str());
					fprintf(fp, "ANALYSIS: %s}\n", rec_id.c_str());
					fprintf(fp, "ANALYSIS: %s}\n", rec_name.c_str());
					fprintf(fp, "ANALYSIS: %s}\n", residue.c_str());
					fprintf(fp, "ANALYSIS: %s}\n", res_id.c_str());
					fprintf(fp, "ANALYSIS: %s}\n\n", chain.c_str());
				}
			}

			fprintf(fp, "DOCKED: MODEL        %d\n", i+1);
			fprintf(fp, "DOCKED: USER    Run = %d\n", i+1);
			fprintf(fp, "DOCKED: USER\n");

			fprintf(fp, "DOCKED: USER    Estimated Free Energy of Binding    =");
			PRINT1000(fp, ((float) (myresults[i].interE + myresults[i].interflexE + torsional_energy)));
			fprintf(fp, " kcal/mol  [=(1)+(2)+(3)-(4)]\n");

			fprintf(fp, "DOCKED: USER\n");

			fprintf(fp, "DOCKED: USER    (1) Final Intermolecular Energy     =");
			PRINT1000(fp, ((float) (myresults[i].interE + myresults[i].interflexE)));
			fprintf(fp, " kcal/mol\n");

			fprintf(fp, "DOCKED: USER        vdW + Hbond + desolv Energy     =");
			PRINT1000(fp, ((float) (myresults[i].interE - myresults[i].interE_elec)));
			fprintf(fp, " kcal/mol\n");

			fprintf(fp, "DOCKED: USER        Electrostatic Energy            =");
			PRINT1000(fp, ((float) myresults[i].interE_elec));
			fprintf(fp, " kcal/mol\n");

			fprintf(fp, "DOCKED: USER        Moving Ligand-Fixed Receptor    =");
			PRINT1000(fp, ((float) myresults[i].interE));
			fprintf(fp, " kcal/mol\n");

			fprintf(fp, "DOCKED: USER        Moving Ligand-Moving Receptor   =");
			PRINT1000(fp, ((float) myresults[i].interflexE));
			fprintf(fp, " kcal/mol\n");

			fprintf(fp, "DOCKED: USER    (2) Final Total Internal Energy     =");
			PRINT1000(fp, ((float) (myresults[i].intraE + myresults[i].intraflexE)));
			fprintf(fp, " kcal/mol\n");

			fprintf(fp, "DOCKED: USER    (3) Torsional Free Energy           =");
			PRINT1000(fp, ((float) torsional_energy));
			fprintf(fp, " kcal/mol\n");

			fprintf(fp, "DOCKED: USER    (4) Unbound System's Energy         =");
			PRINT1000(fp, ((float) (myresults[i].intraE + myresults[i].intraflexE)));
			fprintf(fp, " kcal/mol\n");

			fprintf(fp, "DOCKED: USER\n");
			if(mypars->xml2dlg || mypars->contact_analysis){
				fprintf(fp, "DOCKED: USER    NEWDPF about 0.0 0.0 0.0\n");
				fprintf(fp, "DOCKED: USER    NEWDPF tran0 %.6f %.6f %.6f\n", myresults[i].genotype[0]*mygrid->spacing, myresults[i].genotype[1]*mygrid->spacing, myresults[i].genotype[2]*mygrid->spacing);
				if(!mypars->xml2dlg){
					double phi = myresults[i].genotype[3]/180.0*PI;
					double theta = myresults[i].genotype[4]/180.0*PI;
					fprintf(fp, "DOCKED: USER    NEWDPF axisangle0 %.8f %.8f %.8f %.6f\n", sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta), myresults[i].genotype[5]);
				} else fprintf(fp, "DOCKED: USER    NEWDPF axisangle0 %.8f %.8f %.8f %.6f\n", myresults[i].genotype[3], myresults[i].genotype[4], myresults[i].genotype[5], myresults[i].genotype[GENOTYPE_LENGTH_IN_GLOBMEM-1]);
				fprintf(fp, "DOCKED: USER    NEWDPF dihe0");
				for(j=0; j<myresults[i].reslig_realcoord.num_of_rotbonds; j++)
					fprintf(fp, " %.6f", myresults[i].genotype[6+j]);
				fprintf(fp, "\n");
			}
			fprintf(fp, "DOCKED: USER\n");

			unsigned int lnr=1;
			if ( mypars->flexresfile!=NULL) {
				if ( strlen(mypars->flexresfile)>0 )
					lnr++;
			}
			
			curr_model = pdbqt_template;
			for(atom_cnt = atom_data.size(); atom_cnt-->0;)
			{
				char* line = lineout;
				line += sprintf(line, "%8.3lf", myresults[i].reslig_realcoord.atom_idxyzq[atom_cnt][1]); // x
				line += sprintf(line, "%8.3lf", myresults[i].reslig_realcoord.atom_idxyzq[atom_cnt][2]); // y
				line += sprintf(line, "%8.3lf", myresults[i].reslig_realcoord.atom_idxyzq[atom_cnt][3]); // z
				line += sprintf(line, "%+6.2lf", copysign(fmin(fabs(myresults[i].peratom_vdw[atom_cnt]),99.99),myresults[i].peratom_vdw[atom_cnt])); // vdw
				line += sprintf(line, "%+6.2lf", copysign(fmin(fabs(myresults[i].peratom_elec[atom_cnt]),99.99),myresults[i].peratom_elec[atom_cnt])); // elec
				line += sprintf(line, "    %+6.3lf ", myresults[i].reslig_realcoord.atom_idxyzq[atom_cnt][4]); // q
				line += sprintf(line, "%-2s\n", myresults[i].reslig_realcoord.atom_types[((int) myresults[i].reslig_realcoord.atom_idxyzq[atom_cnt][0])]); // type
				curr_model.insert(atom_data[atom_cnt],lineout);
			}
			fprintf(fp, "%s", curr_model.c_str());
			fprintf(fp, "DOCKED: TER\n");
			fprintf(fp, "DOCKED: ENDMDL\n");
			fprintf(fp, "________________________________________________________________________________\n\n\n");
		}
	}
	// PERFORM CLUSTERING
	// arranging results according to energy, myresults [0] will be the best one (with lowest energy)
	for (j=0; j<num_of_runs-1; j++)
		for (i=num_of_runs-2; i>=j; i--) // arrange according to sum of inter- and intramolecular energies
			if ((myresults [i]).interE+myresults[i].interflexE /*+ (myresults [i]).intraE*/ > (myresults [i+1]).interE+myresults[i+1].interflexE /*+ (myresults [i+1]).intraE*/)	//mimics the behaviour of AD4 unbound_same_as_bound
			//if ((myresults [i]).interE + (myresults [i]).intraE > (myresults [i+1]).interE + (myresults [i+1]).intraE)
			{
				temp_ligres = myresults [i];
				myresults [i] = myresults [i+1];
				myresults [i+1] = temp_ligres;
			}

	for (i=0; i<num_of_runs; i++)
	{
		(myresults [i]).clus_id = 0; // indicates that it hasn't been put into cluster yet
	}

	// the best result is the center of the first cluster
	(myresults [0]).clus_id = 1;
	(myresults [0]).rmsd_from_cluscent = 0;
	num_of_clusters = 1;

	for (i=1; i<num_of_runs; i++) // for each result
	{
		current_clust_center = 0;
		result_clustered = 0;

		for (j=0; j<i; j++) // results with lower id-s are clustered, look for cluster centers
		{
			if ((myresults [j]).clus_id > current_clust_center) // it is the center of a new cluster
			{
				current_clust_center = (myresults [j]).clus_id;
				temp_rmsd = calc_rmsd(&((myresults [j]).reslig_realcoord), &((myresults [i]).reslig_realcoord), mypars->handle_symmetry); // comparing current result with cluster center
				if (temp_rmsd <= cluster_tolerance) // in this case we put result i to cluster with center j
				{
					(myresults [i]).clus_id = current_clust_center;
					(myresults [i]).rmsd_from_cluscent = temp_rmsd;
					result_clustered = 1;
					break;
				}
			}
		}

		if (result_clustered != 1) // if no suitable cluster was found, this is the center of a new one
		{
			num_of_clusters++;
			(myresults [i]).clus_id = num_of_clusters; // new cluster id
			(myresults [i]).rmsd_from_cluscent = 0;
		}
	}

	for (i=1; i<=num_of_clusters; i++) // printing cluster info to file
	{
		subrank = 0;
		cluster_sizes [i-1] = 0;
		sum_energy [i-1] = 0;
		for (j=0; j<num_of_runs; j++)
			if (myresults [j].clus_id == i)
			{
				subrank++;
				(cluster_sizes [i-1])++;
				sum_energy [i-1] += (myresults [j]).interE + myresults[j].interflexE + /*(myresults [j]).intraE +*/ torsional_energy; // intraE can be commented when unbound_same_as_bound
				(myresults [j]).clus_subrank = subrank;
				if (subrank == 1)
				{
					best_energy [i-1] = (myresults [j]).interE + myresults[j].interflexE + /*(myresults [j]).intraE +*/ torsional_energy; // intraE can be commented when unbound_same_as_bound
					best_energy_runid  [i-1] = (myresults [j]).run_number;
				}
			}
	}

	if(mypars->output_dlg){
		// WRITING CLUSTER INFORMATION
		fprintf(fp, "    CLUSTERING HISTOGRAM\n    ____________________\n\n\n");
		fprintf(fp, "________________________________________________________________________________\n");
		fprintf(fp, "     |           |     |           |     |\n");
		fprintf(fp, "Clus | Lowest    | Run | Mean      | Num | Histogram\n");
		fprintf(fp, "-ter | Binding   |     | Binding   | in  |\n");
		fprintf(fp, "Rank | Energy    |     | Energy    | Clus|    5    10   15   20   25   30   35\n");
		fprintf(fp, "_____|___________|_____|___________|_____|____:____|____:____|____:____|____:___\n");

		for (i=0; i<num_of_clusters; i++)
		{
			fprintf(fp, "%4d |", i+1);

			if (best_energy[i] > 999999.99)
				fprintf(fp, "%+10.2e", best_energy[i]);
			else
				fprintf(fp, "%+10.2f", best_energy[i]);
			fprintf(fp, " |%4d |", best_energy_runid[i]);

			if (sum_energy[i]/cluster_sizes[i] > 999999.99)
				fprintf(fp, "%+10.2e |", sum_energy[i]/cluster_sizes[i]);
			else
				fprintf(fp, "%+10.2f |", sum_energy[i]/cluster_sizes[i]);

			fprintf(fp, "%4d |", cluster_sizes [i]);

			for (j=0; j<cluster_sizes [i]; j++)
				fprintf(fp, "#");

			fprintf(fp, "\n");
		}

		fprintf(fp, "_____|___________|_____|___________|_____|______________________________________\n\n\n");

		// writing RMSD table

		fprintf(fp, "    RMSD TABLE\n");
		fprintf(fp, "    __________\n\n\n");

		fprintf(fp, "_____________________________________________________________________\n");
		fprintf(fp, "     |      |      |           |         |                 |\n");
		fprintf(fp, "Rank | Sub- | Run  | Binding   | Cluster | Reference       | Grep\n");
		fprintf(fp, "     | Rank |      | Energy    | RMSD    | RMSD            | Pattern\n");
		fprintf(fp, "_____|______|______|___________|_________|_________________|___________\n" );

		for (i=0; i<num_of_clusters; i++) // printing cluster info to file
		{
			for (j=0; j<num_of_runs; j++)
				if (myresults [j].clus_id == i+1) {
					if (myresults[j].interE + myresults[j].interflexE + torsional_energy > 999999.99)
						fprintf(fp, "%4d   %4d   %4d  %+10.2e  %8.2f  %8.2f           RANKING\n",
						             (myresults [j]).clus_id,
						                   (myresults [j]).clus_subrank,
						                         (myresults [j]).run_number,
						                              myresults[j].interE + myresults[j].interflexE + torsional_energy,
						                                       (myresults [j]).rmsd_from_cluscent,
						                                              (myresults [j]).rmsd_from_ref);
					else
						fprintf(fp, "%4d   %4d   %4d  %+10.2f  %8.2f  %8.2f           RANKING\n",
						             (myresults [j]).clus_id,
						                   (myresults [j]).clus_subrank,
						                         (myresults [j]).run_number,
						                              myresults[j].interE + myresults[j].interflexE + torsional_energy,
						                                       (myresults [j]).rmsd_from_cluscent,
						                                              (myresults [j]).rmsd_from_ref);
				}
		}

		// Add execution and idle time information
		fprintf(fp, "\nRun time %.3f sec", exec_time);
		fprintf(fp, "\nIdle time %.3f sec\n", idle_time);

		if(!mypars->dlg2stdout){
			fclose(fp);
		}
	}

	// if xml has to be generated
	if (mypars->output_xml)
	{
		char* xml_file_name = (char*)malloc(len*sizeof(char));
		strcpy(xml_file_name, mypars->resname);
		strcat(xml_file_name, ".xml");
		fp_xml = fopen(xml_file_name, "w");

		fprintf(fp_xml, "<?xml version=\"1.0\" ?>\n");
		fprintf(fp_xml, "<autodock_gpu>\n");
		fprintf(fp_xml, "\t<version>%s</version>\n",VERSION);
		if((*argc)>1){
			fprintf(fp_xml, "\t<arguments>");
			for(i=1; i<(*argc); i++)
				fprintf(fp_xml, "%s%s", (i>1)?" ":"", argv[i]);
			fprintf(fp_xml, "</arguments>\n");
		}
		if(mypars->dpffile)
			fprintf(fp_xml, "\t<dpf>%s</dpf>\n",mypars->dpffile);
		if(mypars->list_nr>1)
			fprintf(fp_xml, "\t<list_nr>%u</list_nr>\n",mypars->list_nr);
		fprintf(fp_xml, "\t<grid>%s</grid>\n", mypars->fldfile);
		fprintf(fp_xml, "\t<ligand>%s</ligand>\n", mypars->ligandfile);
		if(mypars->flexresfile)
			fprintf(fp_xml, "\t<flexres>%s</flexres>\n",mypars->flexresfile);
		fprintf(fp_xml, "\t<seed>");
		if(!mypars->seed[2]){
			if(!mypars->seed[1]){
				fprintf(fp_xml,"%d", mypars->seed[0]);
			} else fprintf(fp_xml,"%d %d", mypars->seed[0], mypars->seed[1]);
		} else fprintf(fp_xml,"%d %d %d", mypars->seed[0], mypars->seed[1], mypars->seed[2]);
		fprintf(fp_xml, "</seed>\n");
		fprintf(fp_xml, "\t<ls_method>%s</ls_method>\n",mypars->ls_method);
		fprintf(fp_xml, "\t<autostop>%s</autostop>\n",mypars->autostop ? "yes" : "no");
		fprintf(fp_xml, "\t<heuristics>%s</heuristics>\n",mypars->use_heuristics ? "yes" : "no");
		fprintf(fp_xml, "\t<run_requested>%lu</run_requested>\n",mypars->num_of_runs);
		fprintf(fp_xml, "\t<runs>\n");
		double phi, theta;
		for(j=0; j<num_of_runs; j++){
			fprintf(fp_xml, "\t\t<run id=\"%d\">\n",(myresults [j]).run_number);
			if(mypars->contact_analysis){
				if(myresults[j].analysis.size()>0){
					// sort by analysis type
					AnalysisData temp;
					for(unsigned int i=0; i<myresults[j].analysis.size(); i++)
						for(unsigned int k=0; k<myresults[j].analysis.size()-i-1; k++)
							if(myresults[j].analysis[k].type>myresults[j].analysis[k+1].type){
								temp = myresults[j].analysis[k];
								myresults[j].analysis[k]   = myresults[j].analysis[k+1];
								myresults[j].analysis[k+1] = temp;
							}
					fprintf(fp_xml, "\t\t\t<contact_analysis count=\"%lu\">\n", myresults[j].analysis.size());
					std::string types;
					std::string lig_id;
					std::string ligname;
					std::string rec_id;
					std::string rec_name;
					std::string residue;
					std::string res_id;
					std::string chain;
					char item[8], pad[8];
					for(unsigned int i=0; i<myresults[j].analysis.size(); i++){
						if(i>0){
							types    += ",";
							lig_id   += ",";
							ligname  += ",";
							rec_id   += ",";
							rec_name += ",";
							residue  += ",";
							res_id   += ",";
							chain    += ",";
						}
						switch(myresults[j].analysis[i].type){
							case 0: types += "   \"R\"";
							        break;
							case 1: types += "   \"H\"";
							        break;
							default:
							case 2: types += "   \"V\"";
							        break;
						}
						sprintf(item, "%5d ", myresults[j].analysis[i].lig_id);   lig_id+=item;
						sprintf(item, "\"%s\"", myresults[j].analysis[i].lig_name); sprintf(pad, "%6s", item); ligname+=pad;
						sprintf(item, "%5d ", myresults[j].analysis[i].rec_id);   rec_id+=item;
						sprintf(item, "\"%s\"", myresults[j].analysis[i].rec_name); sprintf(pad, "%6s", item); rec_name+=pad;
						sprintf(item, "\"%s\"", myresults[j].analysis[i].residue); sprintf(pad, "%6s", item);  residue+=pad;
						sprintf(item, "%5d ", myresults[j].analysis[i].res_id);   res_id+=item;
						sprintf(item, "\"%s\"", myresults[j].analysis[i].chain); sprintf(pad, "%6s", item);    chain+=pad;
					}
					fprintf(fp_xml, "\t\t\t\t<contact_analysis_types>  %s</contact_analysis_types>\n", types.c_str());
					fprintf(fp_xml, "\t\t\t\t<contact_analysis_ligid>  %s</contact_analysis_ligid>\n", lig_id.c_str());
					fprintf(fp_xml, "\t\t\t\t<contact_analysis_ligname>%s</contact_analsyis_ligname>\n", ligname.c_str());
					fprintf(fp_xml, "\t\t\t\t<contact_analysis_recid>  %s</contact_analysis_recid>\n", rec_id.c_str());
					fprintf(fp_xml, "\t\t\t\t<contact_analysis_recname>%s</contact_analysis_recname>\n", rec_name.c_str());
					fprintf(fp_xml, "\t\t\t\t<contact_analysis_residue>%s</contact_analysis_residue>\n", residue.c_str());
					fprintf(fp_xml, "\t\t\t\t<contact_analysis_resid>  %s</contact_analysis_resid>\n", res_id.c_str());
					fprintf(fp_xml, "\t\t\t\t<contact_analysis_chain>  %s</contact_analysis_chain>\n", chain.c_str());
					fprintf(fp_xml, "\t\t\t</contact_analysis>\n");
				}
			}
			fprintf(fp_xml, "\t\t\t<free_NRG_binding>   %.2f</free_NRG_binding>\n", myresults[j].interE + myresults[j].interflexE + torsional_energy);
			fprintf(fp_xml, "\t\t\t<final_intermol_NRG> %.2f</final_intermol_NRG>\n", myresults[j].interE + myresults[j].interflexE);
			fprintf(fp_xml, "\t\t\t<internal_ligand_NRG>%.2f</internal_ligand_NRG>\n", myresults[j].intraE + myresults[j].intraflexE);
			fprintf(fp_xml, "\t\t\t<torsonial_free_NRG> %.2f</torsonial_free_NRG>\n", torsional_energy);
			fprintf(fp_xml, "\t\t\t<tran0>%.6f %.6f %.6f</tran0>\n", myresults[j].genotype[0]*mygrid->spacing, myresults[j].genotype[1]*mygrid->spacing, myresults[j].genotype[2]*mygrid->spacing);
			phi = myresults[j].genotype[3]/180.0*PI;
			theta = myresults[j].genotype[4]/180.0*PI;
			fprintf(fp_xml, "\t\t\t<axisangle0>%.8f %.8f %.8f %.6f</axisangle0>\n", sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta), myresults[j].genotype[5]);
			fprintf(fp_xml, "\t\t\t<ndihe>%d</ndihe>\n", myresults[j].reslig_realcoord.num_of_rotbonds);
			fprintf(fp_xml, "\t\t\t<dihe0>");
			for(i=0; i<myresults[j].reslig_realcoord.num_of_rotbonds; i++)
				fprintf(fp_xml, "%s%.6f", (i>0)?" ":"", myresults[j].genotype[6+i]);
			fprintf(fp_xml, "\n\t\t\t</dihe0>\n");
			fprintf(fp_xml, "\t\t</run>\n");
		}
		fprintf(fp_xml, "\t</runs>\n");
		fprintf(fp_xml, "\t<result>\n");
		
		fprintf(fp_xml, "\t\t<clustering_histogram>\n");
		for (i=0; i<num_of_clusters; i++)
		{
			fprintf(fp_xml, "\t\t\t<cluster cluster_rank=\"%d\" lowest_binding_energy=\"%.2lf\" run=\"%d\" mean_binding_energy=\"%.2lf\" num_in_clus=\"%d\" />\n",
					i+1, best_energy[i], best_energy_runid[i], sum_energy[i]/cluster_sizes[i], cluster_sizes [i]);
		}
		fprintf(fp_xml, "\t\t</clustering_histogram>\n");
		
		fprintf(fp_xml, "\t\t<rmsd_table>\n");
		for (i=0; i<num_of_clusters; i++)
		{
			for (j=0; j<num_of_runs; j++)
				if (myresults [j].clus_id == i+1)
				{
					fprintf(fp_xml, "\t\t\t<run rank=\"%d\" sub_rank=\"%d\" run=\"%d\" binding_energy=\"%.2lf\" cluster_rmsd=\"%.2lf\" reference_rmsd=\"%.2lf\" />\n",
					                     (myresults [j]).clus_id, (myresults [j]).clus_subrank, (myresults [j]).run_number, myresults[j].interE + myresults[j].interflexE + torsional_energy, (myresults [j]).rmsd_from_cluscent, (myresults [j]).rmsd_from_ref);
				}
		}
		fprintf(fp_xml, "\t\t</rmsd_table>\n");
		fprintf(fp_xml, "\t</result>\n");
		fprintf(fp_xml, "</autodock_gpu>\n");
		fclose(fp_xml);
		free(xml_file_name);
	}
}

void process_result(
                    const Gridinfo*        mygrid,
                    const float*           cpu_floatgrids,
                    const Dockpars*        mypars,
                    const Liganddata*      myligand_init,
                    const Liganddata*      myxrayligand,
                    const int*             argc,
                          char**           argv,
                          SimulationState& sim_state
                   )
{
	std::vector<Ligandresult> cpu_result_ligands(mypars->num_of_runs);

	// Fill in cpu_result_ligands
	float best_energy_of_all = 1000000000000.0;
	for (unsigned long run_cnt=0; run_cnt < mypars->num_of_runs; run_cnt++)
	{
		arrange_result(sim_state.cpu_populations.data()+run_cnt*mypars->pop_size*GENOTYPE_LENGTH_IN_GLOBMEM, sim_state.cpu_energies.data()+run_cnt*mypars->pop_size, mypars->pop_size);
		make_resfiles(sim_state.cpu_populations.data()+run_cnt*mypars->pop_size*GENOTYPE_LENGTH_IN_GLOBMEM,
		              sim_state.cpu_energies.data()+run_cnt*mypars->pop_size,
		              &(sim_state.myligand_reference),
		              myligand_init,
		              myxrayligand,
		              mypars,
		              sim_state.cpu_evals_of_runs[run_cnt],
		              sim_state.generation_cnt,
		              mygrid,
		              cpu_floatgrids,
		              argc,
		              argv,
		              /*1*/0,
		              run_cnt,
		              best_energy_of_all,
		              &(cpu_result_ligands [run_cnt]));
	}

	// Do clustering analysis and generate dlg file
	clusanal_gendlg(cpu_result_ligands.data(),
	                mypars->num_of_runs,
	                myligand_init,
	                mypars,
	                mygrid,
	                argc,
	                argv,
	                sim_state.sec_per_run,
	                sim_state.generation_cnt,
	                sim_state.total_evals/mypars->num_of_runs,
	                sim_state.exec_time,
	                sim_state.idle_time);
}
