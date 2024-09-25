/*

 $Id$

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

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

#include "conformation_sampler.h"
#include "hybrids.h"
#include "ranlib.h"
#include "rep_constants.h"
#include <math.h>

//#define VERBOSE
#define AVOGADRO 6.022e23f
#define RK_CONSTANT 0.0019872065 // constant in entropy calculation
#define TEMP 298 // temperature in entropy calculation
#define RT_CONSTANT 0.592
#define ANGSTROMS_TO_LITERS 1.0e-27f // (1 Angstrom^3 = 10^-27 liters)
#define Vconf 1.0

#define RMSD_SYMMETRY TRUE
#define RMSD_UNIQUE_PAIR TRUE
#define RMSD_HEAVY_ATOMS_ONLY FALSE
#define TRAN_STEP 0.03 // size of translation steps (x,y,z)
#define ROT_ANG_STEP 0.025 // size of step for rotation angle
#define TOR_ANG_STEP 0.03 // size of step for torsion angles
#define RHO 0.5 // parameter for random sampling

#define DEFAULT_RANDOM_SAMPLES 10000
#define DEFAULT_INCREMENTAL_STEPS 5 // note that this is steps up and down, i.e. 
									// X +/- 3

#define ICO_X 0.525731112119133606
#define ICO_Y 0.850650808352039932

//extern class Eval evaluate;

Real (*vt)[SPACE], (*crdpdb)[SPACE];
int (*tlist)[MAX_ATOMS];
Real *lig_center;
int natom;
int *type;
GridMapSetInfo *info;
int h_index;

Real crd[MAX_ATOMS][SPACE];
Real ref_crd[MAX_ATOMS][SPACE];

const Real vertices[12][3] = {{-ICO_X, 0., ICO_Y}, {ICO_X, 0., ICO_Y}, {-ICO_X, 0., -ICO_Y}, {ICO_X, 0., -ICO_Y},
                              {0., ICO_Y, ICO_X}, {0., ICO_Y, -ICO_X}, {0., -ICO_Y, ICO_X}, {0., -ICO_Y, -ICO_X},
                              {ICO_Y, ICO_X, 0.}, {-ICO_Y, ICO_X, 0.}, {ICO_Y, -ICO_X, 0.}, {-ICO_Y, -ICO_X, 0.}};

ConformationSampler::ConformationSampler(const State& init_state, 
 int true_ligand_atoms, Eval *evaluate, int outlev, FILE *logFile) {
	base_state = init_state;
	base_ind = set_ind(info, init_state, evaluate, outlev, logFile);
	base_point = base_ind.phenotyp;	
	base_energy = base_point.evaluate(Normal_Eval);
	cnv_state_to_coords(init_state, vt, tlist, init_state.ntor, crdpdb, base_crd, natom,
	 true_ligand_atoms, outlev, logFile);
	
	dimensionality = BASE_DIMENSIONS + init_state.ntor;
	evals = 0;
	favorable_evals = 0;
	total_energy = 0.0;
	total_favorable_energy = 0.0;
	min_energy = base_energy;
	min_energy_rmsd = 0.0;
	Boltzmann_sum = Boltzmann_diff_sum = 0.0;
	
	// set up the temp variables
	probe_state = base_state;
	probe_ind = base_ind;
	probe_point = base_point;

	// store axis/angle representation in an array

    // read the rep
    Quat base_q = probe_point.readQuat();

	// initialize bounds
	for (int i=0; i < 3; i++) {
		max_values[i] = min_values[i] = base_point.gread(i).real;
		max_values[i+3] = min_values[i+3] = 0.0; // Eulerian angles set to 0
	}
	
	// reset bins
	for (int i=0; i < NUM_BINS; i++) {
		bin_total_energy[i] = 0.0;
		bin_total_favorable_energy[i] = 0.0;
		bin_count[i] = 0;
		bin_count_favorable[i] = 0;
		bin_min_energy[i] = 0.0;
		bin_max_energy[i] = base_energy;
		bin_Boltzmann_sum[i] = 0.0;
	}
}

void ConformationSampler::random_sample(int true_ligand_atoms, int outlev, FILE *logFile) {
	random_sample(1, true_ligand_atoms, outlev, logFile);
}

void ConformationSampler::random_sample(const int num_samples, int true_ligand_atoms, int outlev, FILE *logFile) {
	Real multiplier;
	
	for (int sample=0; sample < num_samples; sample++) {
		probe_point = base_point;
		//multiplier = ranf();
		multiplier = genunf(0.0, 1.25);

		// perturb translation and torsion angles randomly
		for (unsigned int i=0; i < (unsigned int) dimensionality; i++) {
			if (is_rotation_index(i)) continue;
			probe_point.write(probe_point.gread(i).real + gennor(0.0, multiplier*RHO) , i);
			//probe_point.write(probe_point.gread(i).real + genunf(-1.0 * RHO, RHO) , i);
		}

		Real random_axis_angle[4];
		Real new_axis_angle[4];
		const double angle = gennor(0.0, PI/24.0);
		rand_axis(random_axis_angle, angle);
		multiplyraa(base_axis_angle, random_axis_angle, new_axis_angle);

        probe_point.writeQuat(raaToQuat(  new_axis_angle, new_axis_angle[3]));

		current_energy(true_ligand_atoms,outlev,logFile);
	}
}

ConformationSampler::~ConformationSampler(void) {
}

// NOTE: currently, the torsional free energy penalty is not included.
// Since this is an entropic term, I believe we can ignore it in this analysis.
Real ConformationSampler::current_energy(int true_ligand_atoms, int outlev, FILE *logFile) /* not const */ {
	evals++;
	Real energy = probe_point.evaluate(Normal_Eval);
	Real rmsd = current_rmsd(true_ligand_atoms, outlev, logFile);
	
#ifdef VERBOSE
		fprintf(logFile, "state %d %.3f %.3f", evals, energy, rmsd);
		for (int i=0; i < dimensionality; i++) {
			fprintf(logFile, " %.3f", probe_point.gread(i).real);
		}
		fprintf(logFile,"\n");
#endif
	
	total_energy += energy;
	Boltzmann_sum += exp(-energy/RT_CONSTANT);
	Boltzmann_diff_sum += exp(-1.0*(energy - base_energy)/RT_CONSTANT);
	
	// store information on minimum energy conformation
	if (energy < min_energy) {
		min_energy = energy;
		min_energy_rmsd = rmsd;
	}
	
	int bin = (int)(rmsd/BIN_SIZE);
	if (bin < NUM_BINS) {
		bin_count[bin]++;
		bin_total_energy[bin] += energy;
		bin_Boltzmann_sum[bin] += exp(-energy/RT_CONSTANT);
		if (energy < bin_min_energy[bin]) bin_min_energy[bin] = energy;
		if (energy > bin_max_energy[bin]) bin_max_energy[bin] = energy;
	}
	
	if (energy < 0.0) {
		favorable_evals++;
		total_favorable_energy += energy;
		update_bounds();
		
		if (bin < NUM_BINS) {
			bin_count_favorable[bin]++;
			bin_total_favorable_energy[bin] += energy;
		}
	}
	return energy;
}

Real ConformationSampler::current_rmsd(
 int true_ligand_atoms, int outlev, FILE *logFile) /* not const */ {
	probe_ind.phenotyp = probe_point;
	probe_ind.inverse_mapping();
	probe_state = probe_ind.state(base_state.ntor);
	cnv_state_to_coords(probe_state, vt, tlist, probe_state.ntor, crdpdb, crd, natom,
	 true_ligand_atoms, outlev, logFile);
	return getrms(crd, base_crd, RMSD_SYMMETRY, RMSD_UNIQUE_PAIR, natom, type, RMSD_HEAVY_ATOMS_ONLY, h_index);
}

Real ConformationSampler::reference_rmsd(void) const {
	return getrms(base_crd, ref_crd, RMSD_SYMMETRY, RMSD_UNIQUE_PAIR, natom, type, RMSD_HEAVY_ATOMS_ONLY, h_index);
}

void ConformationSampler::update_bounds(void) /* not const */ {
	Real euler[3];
	Real raa[4];
	Real current_value;
    Quat q;
	
    // read the quaternion
    q = probe_point.readQuat();

    // convert to axis-angle
    AxisAngle aa = QuatToAxisAngle( q );

	// set up axis-angle array
    raa[0] = aa.nx;
    raa[1] = aa.ny;
    raa[2] = aa.nz;
    raa[3] = aa.ang;

	raaEuler(raa, euler);
	
	// check existing translation bounds
	for (int i=0; i < 3; i++) {
		current_value = probe_point.gread(i).real;
		if (current_value < min_values[i]) min_values[i] = current_value;
		else if (current_value > max_values[i]) max_values[i] = current_value; 
	}
	
	// check rotation bounds
	for (int i=0; i < 3; i++) {
		if (euler[i] < min_values[i+3]) min_values[i+3] = euler[i];
		else if (euler[i] > max_values[i+3]) max_values[i+3] = euler[i];
	}
}

void ConformationSampler::systematic_search(const int index, int true_ligand_atoms, int outlev, FILE *logFile) {
	
	// for rotation axes, rotate using the pre-defined vertices 
	if ( is_axis_index( index ) ) {
        Quat increment_q;
        Quat new_q;
		
		for (int i=0; i < 12; i++) {
			// set up rotation
            // temp_rotation_angle must have been set before doing this... Potential BUG!
            increment_q = raaToQuat( vertices[i], temp_rotation_angle );
            qmultiply( &new_q, &base_q, &increment_q );
			probe_point.writeQuat( new_q );
			systematic_search(Z_TRANSLATION_INDEX, true_ligand_atoms, outlev, logFile); // go to translation
		}
	}
	
	// translation, rotation angle, and torsion angles
	// step through 
	else {
		int num_steps = DEFAULT_INCREMENTAL_STEPS; // steps up or down
		Real start, step_size;
			
		// set step sizes for different dimensions
		if ( is_translation_index( index ) ) step_size = TRAN_STEP;
		else if ( is_angle_index( index ) ) step_size = ROT_ANG_STEP;
		else step_size = TOR_ANG_STEP;
		
		// for the rotation angle, use different bounds in order to avoid
		// symmetry problems
		if ( is_angle_index( index ) )  start = -2*num_steps*step_size;
		else start = base_point.gread(index).real - num_steps * step_size;
		
		for (int current = 0; current <= 2 * num_steps; current++) {

			if ( is_angle_index( index ) ) {
				temp_rotation_angle = start + current*step_size;
			}
			else {
				probe_point.write(start + current * step_size, index);
			}
			
			if (index == 0) {
                // End recursion
				(void)current_energy(true_ligand_atoms,outlev,logFile);
			}
			else {
				// check if the rotation angle is 0, to avoid shifting axis unnecessarily
				if ( is_angle_index( index ) && current == 2 * num_steps) {
                    probe_point.writeQuat( base_point.readQuat() );
					systematic_search(Z_TRANSLATION_INDEX, true_ligand_atoms, outlev, logFile); // skip to translation
					//current_energy(true_ligand_atoms,outlev,logFile);// DEBUGGING
				}
				else {
					systematic_search(index-1, true_ligand_atoms, outlev, logFile);
				}
			}
		}
	}
}

Real ConformationSampler::fraction_favorable(void) const {
	return favorable_evals/evals;
}

Real ConformationSampler::average_favorable_energy(void) const {
	if (favorable_evals == 0) return 0;
	else return total_favorable_energy/favorable_evals;
}

Real ConformationSampler::energy_volume(void) const { 
	return total_favorable_energy/evals;
}

Real ConformationSampler::configurational_integral(void) const {
	Real Vb = 1.0;
	for (int i=0; i < 6; i++) {
		Vb *= (max_values[i]-min_values[i]);
	}

	Vb *= ANGSTROMS_TO_LITERS;
	return Vb;
}

/*
 * estimate entropy, as described by Ruvinsky and Kozintsev
 * return (0.0019872065)*(298)*math.log(Vb * 6.02 * 10**23 / (8*math.pi*math.pi))
 */

Real ConformationSampler::RK_entropy(void) const {
	return RK_CONSTANT * TEMP * log(configurational_integral() * AVOGADRO/ (8 * PI * PI));
}

Real ConformationSampler::partition_function(void) const {
	return -RT_CONSTANT*log(Boltzmann_sum/evals);
}

Real ConformationSampler::partition_function(const int bin) const {
	return -RT_CONSTANT*log(bin_Boltzmann_sum[bin]/bin_count[bin]);
}

Real ConformationSampler::normalized_volume(void) const {
	Real volume = 0.0;
	for (int i=0; i < NUM_BINS; i++) {
		volume += bin_total_favorable_energy[i]/bin_count[i]/NUM_BINS;
	}
	return volume;
}

Real ConformationSampler::normalized_Boltzmann(void) const {
	Real boltzmann_sum = 0.0;
	for (int i=0; i < NUM_BINS; i++) {
		boltzmann_sum += partition_function(i)/NUM_BINS;
	}
	return boltzmann_sum;
}

Real ConformationSampler::entropy_estimate(void) const {
	Real Vtot = AVOGADRO/(8*PI*PI);
	Vtot *= pow(1/(2*PI), (dimensionality - BASE_DIMENSIONS));
	//fprintf(logFile, "Vtot: %g\n", Vtot);
	return RT_CONSTANT*log(Vtot*Vconf*Boltzmann_diff_sum);
}

void ConformationSampler::output_statistics(int outlev, FILE *logFile) const {
	if(outlev<0) return;
	fprintf(logFile, "Conformation starting energy: %.3f\n", base_energy);
	fprintf(logFile, "RMSD from reference state: %.3f\n", reference_rmsd());
	fprintf(logFile, "Fraction of favorable evaluations: %.3f\n", (Real)favorable_evals/evals);
	fprintf(logFile, "Average favorable energy: %.3f\n", total_favorable_energy/favorable_evals);
	fprintf(logFile, "Estimated energy volume: %.3f\n", total_favorable_energy/evals);
	//fprintf(logFile, "Normalized estimated energy volume: %.3f\n", normalized_volume());
	fprintf(logFile, "Vb estimate: %.3f\n", Vconf*Boltzmann_diff_sum);
	fprintf(logFile, "Entropy estimate: %.3f\n", entropy_estimate());
	fprintf(logFile, "Boltzmann-weighted energy: %.3f\n", partition_function());
	//fprintf(logFile, "Normalized Boltzmann-weighted energy: %.3f\n", normalized_Boltzmann());
	fprintf(logFile, "Minimum energy found: %.3f (%.3f A from starting point)\n", min_energy, min_energy_rmsd);
	//fprintf(logFile, "Bins in local region.\n");
	fprintf(logFile, "\nRMSD       #     fraction  Volume    Avg. (-)   Min E     Max E    Boltzmann\n");
	for (int i=0; i < NUM_BINS; i++) {
		fprintf(logFile, "%.1f    %7d    %2.3f    %2.3f    %2.3f    %2.3f    %2.3f    %2.3f\n", (i+1)*BIN_SIZE, bin_count[i], (Real)bin_count_favorable[i]/bin_count[i], bin_total_favorable_energy[i]/bin_count[i], bin_total_favorable_energy[i]/bin_count_favorable[i], bin_min_energy[i], bin_max_energy[i], partition_function(i));
	}
	fprintf(logFile, "%d evaluations.\n\n", evals);
}

void systematic_conformation_sampler(const State hist[MAX_RUNS], const int nconf, Real init_vt[MAX_TORS][SPACE], Real init_crdpdb[MAX_ATOMS][SPACE], int init_tlist[MAX_TORS+1][MAX_ATOMS], Real init_lig_center[SPACE], const int init_natom, int init_type[MAX_ATOMS], GridMapSetInfo *const init_info,
 int true_ligand_atoms, Eval *evaluate, int outlev, FILE *logFile) {
	vt = init_vt;
	crdpdb = init_crdpdb;
	tlist = init_tlist;
	lig_center = init_lig_center;
	natom = init_natom;
	type = init_type;
	info = init_info;
	
	setup_reference_coordinates();
	
	fprintf(logFile, "Initiating a systematic search.\n");
	for (int i=0; i < nconf; i++) {
		fprintf(logFile, "\nConformation %d:\n", i+1);
		State base_state = hist[i];
		ConformationSampler CS(base_state, true_ligand_atoms, evaluate, outlev, logFile);
		//CS.systematic_search(CS.dimensionality-1, true_ligand_atoms, outlev, logFile);
		CS.systematic_search(BASE_DIMENSIONS-1, true_ligand_atoms, outlev, logFile);
		CS.output_statistics(outlev, logFile);
	}
	fprintf(logFile,"\n\n");
}

void random_conformation_sampler(const State hist[MAX_RUNS], const int nconf, /* not const */ int num_samples, Real init_vt[MAX_TORS][SPACE], Real init_crdpdb[MAX_ATOMS][SPACE], int init_tlist[MAX_TORS+1][MAX_ATOMS], Real init_lig_center[SPACE], const int init_natom, int init_type[MAX_ATOMS], GridMapSetInfo *const init_info,
 int true_ligand_atoms, Eval *evaluate, int outlev, FILE *logFile) {
	vt = init_vt;
	crdpdb = init_crdpdb;
	tlist = init_tlist;
	lig_center = init_lig_center;
	natom = init_natom;
	type = init_type;
	info = init_info;
	
	setup_reference_coordinates();
	
	if (num_samples == 0) num_samples = DEFAULT_RANDOM_SAMPLES;
	
	fprintf(logFile, "Initiating a random search using %d samples near each conformation.\n", num_samples);
	for (int i=0; i < nconf; i++) {
		fprintf(logFile, "\nConformation %d:\n", i+1);
		State base_state = hist[i];
		ConformationSampler CS(base_state, true_ligand_atoms, evaluate, outlev, logFile);
		CS.random_sample(num_samples, true_ligand_atoms, outlev, logFile);
		CS.output_statistics(outlev, logFile);
	}
	
	fprintf(logFile,"\n\n");
}


/* copied (and slightly modified) from non-included code in call_glss.cc */
Individual set_ind(GridMapSetInfo *const info, const State state, Eval *evaluate, int outlev, FILE *logFile)
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;
   int i;

   temp_Gtype = generate_Gtype(state.ntor, info, outlev, logFile);
   temp_Ptype = generate_Ptype(state.ntor, info, evaluate, outlev, logFile);

   // use the state to generate a Genotype
   temp_Gtype.write( state.T.x, 0 );
   temp_Gtype.write( state.T.y, 1 );
   temp_Gtype.write( state.T.z, 2 );
   temp_Gtype.write( state.Q.x, 3 );
   temp_Gtype.write( state.Q.y, 4 );
   temp_Gtype.write( state.Q.z, 5 );
   temp_Gtype.write( state.Q.w, 6 );
   for (i=0;i<state.ntor; i++) {
       temp_Gtype.write( state.tor[i], 7+i );
   };

   Individual temp(temp_Gtype, temp_Ptype);   

   // use mapping to generate a Phenotype
   //temp.phenotyp = temp.mapping();
   temp.mapping();

   return(temp);
}

void raaEuler(const Real raa[4], /* not const */ Real euler[3]) {
	Real s = sin(raa[3]);
	Real c = cos(raa[3]);
	Real t = 1.0 - c;
	
	// check for singularities
	if (raa[0]*raa[1]*t + raa[2]*s > 0.998) {
		euler[0] = 0.0;
		euler[1] = atan2(raa[0]*sin(raa[3]/2), cos(raa[3]/2));
		euler[2] = PI/2;
	}
	
	else if (raa[0]*raa[1]*t + raa[2]*s < -0.998) {
		euler[0] = 0.0;
		euler[1] = -atan2(raa[0]*sin(raa[3]/2), cos(raa[3]/2));
		euler[2] = -PI/2;
	}
	
	euler[0] = atan2(raa[0]*s - raa[1]*raa[2]*t , 1 - (raa[0]*raa[0] + raa[2]*raa[2])*t);
	euler[1] = atan2(raa[1]*s - raa[0]*raa[2]*t , 1 - (raa[1]*raa[1] + raa[2]*raa[2])*t);
	euler[2] = asin(raa[0]*raa[1]*t + raa[2]*s);
}

void raaMatrix(/* not const */ Real raa[4], /* not const */ Real matrix[3][3]) {
	Real angle_cos = cos(raa[3]);
	Real angle_sin = sin(raa[3]);
	Real t = 1.0 - angle_cos;
	
	// make sure that input vecotr is a unit vector
	Real length = hypotenuse(raa[0], raa[1], raa[2]);
	raa[0] /= length;
	raa[1] /= length;
	raa[2] /= length;
	
	matrix[0][0] = angle_cos + raa[0]*raa[0]*t;
	matrix[1][1] = angle_cos + raa[1]*raa[1]*t;
	matrix[2][2] = angle_cos + raa[2]*raa[2]*t;
	
	Real tmp1 = raa[0]*raa[1]*t;
	Real tmp2 = raa[2]*angle_sin;
	matrix[1][0] = tmp1 + tmp2;
	matrix[0][1] = tmp1 - tmp2;
	
	tmp1 = raa[0]*raa[2]*t;
	tmp2 = raa[1]*angle_sin;
	matrix[2][0] = tmp1 - tmp2;
	matrix[0][2] = tmp1 + tmp2;
	
	tmp1 = raa[1]*raa[2]*t;
	tmp2 = raa[0]*angle_sin;
	matrix[2][1] = tmp1 + tmp2;
	matrix[1][2] = tmp1 - tmp2;
}

void matrixraa(const Real matrix[3][3], /* not const */ Real raa[4]) {
	Real length = hypotenuse(matrix[2][1] - matrix[1][2], matrix[2][0] - matrix[0][2], matrix[1][0] - matrix[0][1]);
	
	// need to check acos() parameter to avoid values out of range
	Real cosine = (matrix[0][0] + matrix[1][1] + matrix[2][2] - 1)/2;
	if (cosine > 1.0) raa[3] = 0;
	else if (cosine < -1.0) raa[3] = PI;
	else raa[3] = acos(cosine);
	
	raa[0] = (matrix[2][1] - matrix[1][2])/length;
	raa[1] = (matrix[0][2] - matrix[2][0])/length;
	raa[2] = (matrix[1][0] - matrix[0][1])/length;
}

void multiplyraa(/* not const */ Real raa1[4], /* not const */ Real raa2[4], /* not const */ Real raa_result[4]) {
	Real matrix1[3][3];
	Real matrix2[3][3];
	Real result_matrix[3][3];
	
	raaMatrix(raa1, matrix1);
	raaMatrix(raa2, matrix2);
	matrixMultiply(matrix1, matrix2, result_matrix);
	matrixraa(result_matrix, raa_result);
}

void matrixMultiply(const Real m1[3][3], const Real m2[3][3], /* not const */ Real result[3][3]) {
	result[0][0] = m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0] + m1[0][2]*m2[2][0];
	result[0][1] = m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1] + m1[0][2]*m2[2][1];
	result[0][2] = m1[0][0]*m2[0][2] + m1[0][1]*m2[1][2] + m1[0][2]*m2[2][2];
	result[1][0] = m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0] + m1[1][2]*m2[2][0];
	result[1][1] = m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1] + m1[1][2]*m2[2][1];
	result[1][2] = m1[1][0]*m2[0][2] + m1[1][1]*m2[1][2] + m1[1][2]*m2[2][2];
	result[2][0] = m1[2][0]*m2[0][0] + m1[2][1]*m2[1][0] + m1[2][2]*m2[2][0];
	result[2][1] = m1[2][0]*m2[0][1] + m1[2][1]*m2[1][1] + m1[2][2]*m2[2][1];
	result[2][2] = m1[2][0]*m2[0][2] + m1[2][1]*m2[1][2] + m1[2][2]*m2[2][2];
}

void rand_axis(/* not const */ Real axis[4], const double angle) {
	axis[2] = genunf(-1.0, 1.0);
	const Real t = genunf(0.0, 2*PI);
	const Real w = sqrt(1 - axis[2]*axis[2]);
	axis[0] = w * cos(t);
	axis[1] = w * sin(t);
	axis[3] = (Real) angle;
}

void setup_reference_coordinates(void) {
	for (int i = 0;  i < natom;  i++) {
		ref_crd[i][0] = lig_center[0] + crdpdb[i][0];
		ref_crd[i][1] = lig_center[1] + crdpdb[i][1];
		ref_crd[i][2] = lig_center[2] + crdpdb[i][2];
	}
}
