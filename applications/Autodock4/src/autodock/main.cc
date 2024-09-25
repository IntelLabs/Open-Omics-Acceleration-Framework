/* AutoDock
 $Id: main.cc,v 1.213 2014/07/10 23:49:14 mp Exp $

**  Function: Performs Automated Docking of Small Molecule into Macromolecule
**Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
** All Rights Reserved.
**____________________________________________________________________________
** Primary Authors: 
**            Garrett Matthew Morris, C/C++ version 
**            David Goodsell, Original FORTRAN version 1.0
**                                       e-mail: goodsell@scripps.edu
**
** Other Contributors: see file AUTHORS
**            
**            Laboratory of Arthur J. Olson
**            The Scripps Research Institute
**            Department of Molecular Biology, MB5
**            10550 North Torrey Pines Road
**            La Jolla, CA 92037.
**____________________________________________________________________________
**    Inputs: Docking parameter file, Small Molecule PDBQT file, 
**  macromolecular grid map files.
**   Returns: Autodock Log File, includes docked conformation clusters (PDBQT)
**____________________________________________________________________________
** Modification Record (pre-CVS)
** Date     Inits   Comments
** 09/06/95 RSH     Added code to handle GA/LS stuff
**          DSG     Quaternion rotations
**          DSG     Generates torsions from annotated pdb list
**          DSG     Generates internal energies
**          DSG     Performs a limited Cluster analysis of conformations
** 05/07/92 GMM     C translation
** 05/14/92 GMM     Time-dependent seed in random-number generation
** 10/29/92 GMM     Application Visualization System (AVS) readable grid
**                  display file input.
**                  [AVS is a trademark of Stardent Computer Inc.]
**                  Also added the 'total_charge' check.
** 11/19/93 GMM     #ifdef NOSQRT, with non-square-rooting acceleration.
** 09/26/94 GMM     Cluster analysis now outputs RMS deviations.
** 09/28/94 GMM     Modularized code.
** 10/02/94 GMM     Distance constraints added, for Ed Moret. Accelerated.
** 09/06/95 RSH     Incorporation of GA/SW tokens
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
*******************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <sys/types.h> // time_t time(time_t *tloc);
#include <time.h>      // time_t time(time_t *tloc);
#include <stdlib.h>
#ifndef HAVE_SYSCONF
#include "mingw_sysconf.h"  // for sysconf(_SC_CLK_TCK) and possibly gethostname
#endif
#include <string>
using std::string;
#define streq(a,b) (0==strcasecmp(a,b)) // case-independent string match

// convenience macro for parsing the (many) single-argument DPF lines:
//  if not 1 (non-ignored) argument, stop, reporting fatal error
#define get1arg(line, fmt, addr, token) if(1!=sscanf(line, fmt, addr))stop("syntax error in " token " line")
// convenience macro for reporting syntax errors in DPF lines:
#define syntaxstop(s) {char ss[LINE_LEN+50];sprintf(ss,"syntax error or illegal value in %s line",s);stop(ss);}

// convenience macro for making value boolean (in place)
#define mkbool(x) (x=((x)!=0))
// convenience macro for plural noun string
#define pl(i) ((i==1)?"":"s")

#include <sys/param.h>
#include <ctype.h> // tolower
#include <unistd.h> // sysconf and getcwd

/* the BOINC API header file */
#ifdef BOINC
#include "diagnostics.h"
#include "boinc_api.h"
#include "filesys.h"                 // boinc_fopen(), etc... */
#endif

#include "coliny.h"
#include "hybrids.h"
#include "ranlib.h"
#include "gs.h"
#include "ls.h"
#include "rep.h"
#include "support.h"
#include "distdepdiel.h"
#include "calculateEnergies.h"
#include "conformation_sampler.h"
#include "main.h"
#include "threadlog.h"
#include "alea.h"
#include "timesys.h" // for struct tms
// PSO
//#include "call_cpso.h"
#include "pso.h"
#include "dimLibrary.h"
#include "center_ligand.h"

/* globals : */
extern int debug;
extern int keepresnum;
extern Real idct;
Eval evaluate; // used by the search methods that are not yet thread-safe 
int sel_prop_count = 0; // gs.cc debug switch


static const char* const ident[] = {ident[1], "@(#)$Id: main.cc,v 1.213 2014/07/10 23:49:14 mp Exp $"};




// static (local to this source file main.cc) DPF-parsing state variables: 

static Boole parameter_library_found = FALSE; // was atom parm file specified? (not required)
static Boole B_found_about_keyword = FALSE; //set false by 'move' true by 'about'
static Boole B_found_tran0_keyword = FALSE; //set false by 'move' true by 'tran0'
static Boole B_found_elecmap = FALSE;
static Boole B_found_desolvmap = FALSE;
static Boole B_atom_types_found = FALSE;
static Boole B_havemap = FALSE;
static Boole B_found_move_keyword = FALSE;
static Boole B_found_ligand_types = FALSE;
static Boole B_found_autodock_parameter_version = FALSE;
static Boole B_have_flexible_residues = FALSE;  // does the receptor have flexible residues

static int true_ligand_atoms = 0; // used by exit_if ... 



 /* local-to-main functions: */
static void exit_if_missing_elecmap_desolvmap_about(string keyword); // see bottom of main.cc
static int getoutlev(char *line, int *outlev); // see bottom of main.cc  0==fail, 1==OK
static void set_seeds( FourByteLong seed[2], char seedIsSet[2], FourByteLong runseed[][2], const int outlev, FILE *logFile ); // see below
static int processid();

// PSO  - Particle Swarm Optimization  - not officially supported
// State Structure Variable DECLARATION
int S ; // Swarm size
double pso_xmin[PSO_D_MAX], pso_xmax[PSO_D_MAX]; // Intervals defining the search space

#ifdef _OPENMP
/* M Pique */
#include <omp.h>
#else
#define omp_get_thread_num() (0)
#define omp_get_max_threads() (1)
#endif


int main (int argc, const char ** argv)


{

//   MAX_GRID_PTS & MAX_MAPS
//
static MapType *map;  // Use this with malloc, see grid.h and map_declare.h
// map is used as map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];

GridMapSetInfo *info;  // this information is from the AVS field file

//   MAX_ATOMS
//
char atomstuff[MAX_ATOMS][MAX_CHARS];
char pdbaname[MAX_ATOMS][5];
Real crdorig[MAX_ATOMS][SPACE];  // original coords
Real crdpdb[MAX_ATOMS][SPACE];  // PDB coordinates, recentered by "about", possibly reoriented
Real charge[MAX_ATOMS];
Real abs_charge[MAX_ATOMS];
Real qsp_abs_charge[MAX_ATOMS];
int type[MAX_ATOMS];
int bond_index[MAX_ATOMS];
int ignore_inter[MAX_ATOMS];
/* the following are modified according to state: */
Real crd[MAX_ATOMS][SPACE];     // current coordinates according to State
EnergyComponent peratomE[MAX_ATOMS];

//   MAX_TORS
// These are all constant for any given ligand and torsion tree combination
int  tlist[MAX_TORS+1][MAX_ATOMS];
Real vt[MAX_TORS][SPACE];
Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2];
unsigned short  US_TorE[MAX_TORS];
Boole B_isTorConstrained[MAX_TORS];
int N_con[MAX_TORS];
unsigned short US_torProfile[MAX_TORS][NTORDIVS];

//   MAX_NONBONDS
//
NonbondParam *nonbondlist = new NonbondParam[MAX_NONBONDS];

//   LINE_LEN
//
char error_message[LINE_LEN+100];
char message[LINE_LEN];
char line[LINE_LEN];
char torfmt[LINE_LEN];
char param[2][LINE_LEN];
char rms_atoms_cmd[LINE_LEN];
char c_mode_str[LINE_LEN];
char confsampler_type[LINE_LEN];
char autodock_parameter_version[LINE_LEN]; //eg 4.1.1

// filename max length is taken from system include file
char FN_parameter_library[PATH_MAX];
char FN_ligand[PATH_MAX];
char FN_flexres[PATH_MAX];
char FN_rms_ref_crds[PATH_MAX];
char FN_clus[PATH_MAX];
char FN_watch[PATH_MAX];
char dummy_FN_ligand[PATH_MAX];
char FN_pop_file[PATH_MAX];
char FN_trj[PATH_MAX];
char FN_current_working_directory[PATH_MAX];

//   MAX_CHARS
char hostnm[MAX_CHARS];

//   MAX_RECORDS
//
char PDBQT_record[MAX_RECORDS][LINE_LEN];

//   SPACE (3)
//
Real lig_center[SPACE];

//   MAX_RUNS
//
Real econf[MAX_RUNS];  // this is the list of energies printed in the histogram in "analysis"
State sHist[MAX_RUNS];  /* qtnHist[MAX_RUNS][QUAT],torHist[MAX_RUNS][MAX_TORS];*/
FourByteLong runseed[MAX_RUNS][2]; /* initial seed for each run, computed ahead of time */

State sUnbound; // State of the unbound ligand's conformation
State sUnbound_ext; // State of the unbound ligand's conformation after extended-conformation search
//// UNCOMMENT if using Step 2 in unbound calculation ---> State sUnbound_ls; // State of the unbound ligand's conformation after a local search
State sUnbound_ad; // State of the unbound ligand's conformation after an AutoDock search

char S_contype[8];

//   MAX_ATOM_TYPES
//
char            *ligand_atom_type_ptrs[MAX_ATOM_TYPES]; /* array of ptrs used to parse input line of atom type names */
ParameterEntry parameterArray[MAX_ATOM_TYPES]; // input  nonbond and desolvation parameters
static ParameterEntry * foundParameter;

// internal DPF-parsing state variables
static char seedIsSet[2]; // starts out empty (NULL chars)

// gaussian torsion dihedral contraints
Boole B_isGaussTorCon = FALSE;
Boole B_constrain_dist;
int iCon=0;
Real F_A;
Real F_Aova;
Real F_tor;
Real F_torPref;
Real F_torHWdth;

Real sqlower;
Real squpper;

// ELECSCALE converts between CGS units and SI units;
// see, e.g. p 254, "Molecular Modeling and Simulation", by Tamar Schlick, Springer.
//
// Units of ELECSCALE are (Kcal/mol ) * (Angstrom / esu^2)
// and this allows us to use distances in  Angstroms and charges in esu...
const Real ELECSCALE = 332.06363;

// const Real ELECSCALE = 83.0159075;   this ELECSCALE (corresponding to eps(r) = 1/4r) gives -7.13 kcal/mol for 1pgp Tests/test_autodock4.py


//Real c=0.0;

// energy evaluation working storage:
Real cA;
Real cB;
Real torsFreeEnergy = 0.0;
Real AD3_FE_coeff_estat   = 1.000; // obsolete option in intelec 


// ligand setup - these must (should be..) reset for each serial ligand (not really supported - MPique)
Real maxrad = -1.0;
Real r2sum=0.0;
static Boole B_haveCharges=FALSE;


// Simulated annealing  DPF-settable parameters:
Boole B_linear_schedule = TRUE; /* TRUE is ADT default, other is geometric  */
	// ^^ MP TODO 2012 BY making TRUE default, there is no way to turn off
Boole B_selectmin = TRUE; // adopt min instead of last state - ADT default TRUE
Real e0max = 0; // minimum energy for simanneal initial state - 0 is ADT default
int MaxRetries = 10000; // maximum number of retries for simanneal ligand init. 10000 is ADT default
Real RT0 = /* 616.0*/ 100.; /* 616.0 was pre-4.2.5 default */
Real RTFac = 0.90; /* 0.95 was pre-4.2.5 default */
int ncycles = 50; /* 50 is ADT default */
int naccmax = 30000; /* 100 was pre-4.2.5 default */
int nrejmax = 30000; /* 100 was pre-4.2.5 default */

 // note: trnStep0, qtwStep0, torStep0 also control 'investigate'
 // but (despite appearances) do not control GA mutations (MP 2012)
Real trnFac = 1.0; /* 1.0 is ADT default: i.e., no reduction in geom sched */
Real trnStep0 = 0.2;  /* 2 was pre-4.2.5 default */
Real trnStepFinal = 0; // no default value, must be set in DPF
Real qtwFac = 1.0; /* 1.0 is ADT default: i.e., no reduction in geom sched */
Real qtwStep0 = DegreesToRadians( 5.0 );  /* 50 was pre-4.2.5 default */
Real qtwStepFinal = 0; // no default value, must be set in DPF
Real torStep0 = DegreesToRadians( 5.0 );  /* 50 was pre-4.2.5 default */
Real torStepFinal = 0; // no default value, must be set in DPF
Real torFac = 1.0; /* 1.0 is ADT default: i.e., no reduction in geom sched */
 // simanneal file-based or real-time monitoring
Boole B_write_trj = FALSE;
Boole B_watch = FALSE;
Boole B_acconly = FALSE;
Boole B_either = FALSE;

// simanneal internal variables, not directly settable in DPF:
Real RTreduc; // RT decrease per cycle if linear_schedule
// Real RJ = 8.31441;     // in J/K/mol, Gas Constant, Atkins Phys.Chem., 2/e
// Real Rcal = 1.9871917; // in cal/K/mol, Gas Constant, RJ/4.184
// Real T0K = 273.15;        // 0 degrees Celsius, in K
Boole B_tempChange = FALSE;
Boole B_torReduc = FALSE;
Boole B_trnReduc = FALSE;
Boole B_qtwReduc = FALSE;
Boole B_CalcTrnRF = FALSE;
Boole B_CalcQtwRF = FALSE;
Boole B_CalcTorRF = FALSE;



Real unbound_internal_FE = 0.0;
Real unbound_ext_internal_FE = 0.0;
Real unbound_ad_internal_FE = 0.0;


//  LS  local search (SW or PSW)
Real rho = 1.0; // for SW
Real lb_rho = 0.01; // for SW
Real *rho_ptr = NULL; // for PSW array of rho
Real *lb_rho_ptr = NULL; // for PSW array of lb_rho
Real psw_trans_scale = 1.0; // 1 angstrom
Real psw_rot_scale = 0.05;  //about 3 degrees, we think
Real psw_tors_scale = 0.1; //about 6 degrees

// energy evaluation and potential scoring function parameters settable in DPF
Unbound_Model ad4_unbound_model = Unbound_Default; //NOT Same_As_Bound so user can specify in dpf
Real scale_1_4 = 0.5;
Real scale_eintermol = 1.0; // scale factor for intermolecular energy term vs intra
Boole B_include_1_4_interactions = FALSE;  // FALSE was the default behaviour in AutoDock versions 1 to 3.
Boole B_use_non_bond_cutoff = TRUE;
Boole B_calcIntElec = TRUE;
Boole B_calcIntElec_saved = FALSE;
Real r_smooth=0.5; // vdw nonbond smoothing range, not radius, Ang - default 0.5 matches AutoGrid recommendations
Real WallEnergy = 1000; /* Energy barrier beyond walls of gridmaps. 1000 is ADT default  */

int xA = 12;
int xB = 6;
Real cA_unbound = 392586.8;  // repulsive
Real cB_unbound = 0.0; // attractive

// set by read_parameter_library:
Linear_FE_Model AD4;


// Distance-dependence in Desolvation Term
const double sigma = 3.6L;
const double qsolpar = 0.01097L;


EnergyBreakdown eb;
initialise_energy_breakdown(&eb, 0, 0);

static Output_pop_stats output_pop_stats;



// initial population for LS, GA, GALS, PSO
Boole B_RandomTran0 = TRUE;
Boole B_RandomQuat0 = TRUE;
Boole B_RandomDihe0 = TRUE;
Boole B_reorient_random = FALSE; // if true, create a new random orientation before docking





// cluster analysis
Boole B_cluster_mode = FALSE; // if TRUE, writes to file named in DPF "cluster"
Boole B_symmetry_flag = TRUE;
Boole B_unique_pair_flag = FALSE;
Boole B_write_all_clusmem = FALSE;
Boole B_ShowTorE = FALSE;
Boole B_rms_atoms_ligand_only = TRUE;  // cluster on the ligand atoms only
Boole B_rms_heavy_atoms_only = FALSE;  // cluster on the ligand heavy atoms only, exclude hydrogens
Real clus_rms_tol = 2.0; // 2.0 is ADT default

// ligand atom types and matching receptor maps
Boole B_charMap = FALSE;

int atm1=0;
int atm2=0;
int a1=0;
int a2=0;
int atomC1;
int atomC2;
int dpf_keyword = -1;
int h_index = -1; //index of hydrogen type if any 
int n_heavy_atoms_in_ligand = 0;
int indcom = 0;

// affinity/inhibition reporting controls: 
Real torsdoffac = 0.3113;
int ligand_is_inhibitor = 1;

int nruns = 50; // for GA/GALS and SIMANNEAL

int natom = 0;

// For energy breakdown of non-bonded interactions:
int     Nnb_array[3] = {0};    // number of nonbonds in the ligand, intermolecular and receptor groups
static GroupEnergy group_energy; // energy components of each of the five groups (intra-ligand, inter, and intra-receptor...)
Boole B_havenbp = FALSE;

int nconf = 0;  // overall count of number of runs so far, must < MAX_RUNS
int nlig = 0;
int nres = 0;
int nmol = 0;
int Nnb = 0;
int ntor = 0;
int ntor_ligand = 0;
int ntorsdof = 0;
int num_maps = 0;
int num_atom_types = 0;
int nval = 0;
int nfields = 0;
int trj_end_cyc = 0;
int trj_begin_cyc = 0;
int trj_freq = 0;
int xA_unbound = 12;
int xB_unbound = 6;
int I_tor;
int I_torBarrier;

// for INVESTIGATE operation
int OutputEveryNTests = 1000;
int NumLocalTests = 10;
int maxTests = 10000;

/* int beg; */
/* int end; */
/* int imol = 0; */
// unsigned int min_evals_unbound =  250000;
unsigned int max_evals_unbound = 1000000;
int saved_sInit_ntor = 0;
int confsampler_samples = 0;

unsigned short US_energy;
unsigned short US_tD;
unsigned short US_torBarrier = TORBARMAX;
unsigned short US_min = TORBARMAX;

register int i = 0;
register int j = 0;
register int k = 0;
//register int m = 0;
register int xyz = 0;

State sInit;            /* Real qtn0[QUAT], tor0[MAX_TORS]; */

Quat q_reorient;

Molecule ligand;        /* ligand */

static Real F_A_from;
static Real F_A_to;
static Real F_lnH;
static Real F_W;
static Real F_hW;
static long clktck = 0;

#ifndef VERSION
static string version_num = "4.2.2";
#else
static string version_num = VERSION;
#endif

struct tms tms_jobStart;
Clock  jobStart;


EnergyTables *ad_energy_tables;  // Holds vdw+Hb, desolvation & dielectric lookup tables
EnergyTables *unbound_energy_tables;  // Use for computing unbound energy & conformation

Statistics map_stats;

//  GA parameters controlled in DPF
static FourByteLong seed[2]; // also used by simanneal & investigate as of 4.2.5 release (default is process id, time)
unsigned int pop_size = 150; // 150 is ADT default
unsigned int num_generations = 0;  //  Don't terminate on the basis of number of generations
unsigned int num_evals = 250000;
unsigned int num_evals_unbound = num_evals;
Selection_Mode s_mode = Proportional;
int elitism = 1; // 1 is ADT default
Real linear_ranking_selection_probability_ratio = 2.0;
Xover_Mode c_mode = TwoPt;  //  can be: OnePt, TwoPt, Uniform or Arithmetic
Real m_rate = 0.02; // 0.02 is ADT default
Real c_rate = 0.80; // 0.80 is ADT default
Real alpha = 0; // I believe is unused in existing GA code - MP 2012
Real beta = 1;  // I believe is unused in existing GA code - MP 2012
Real localsearch_freq = 0.06;  // 0.06 is ADT default
Worst_Mode w_mode = AverageOfN; // note: no keyword to change this in DPF
int window_size = 10;
int low = 0;  // unsure if is used in existing GA code - MP 2012
int high = 100;  // unsure if is used in existing GA code - MP 2012

// local search (PSW) parameters controlled in DPF
unsigned int max_its = 300;
unsigned int max_succ = 4;
unsigned int max_fail = 4;


// internal variables for GA and GA/LS
// For Branch Crossover Mode
int end_of_branch[MAX_TORS];


EvalMode e_mode = Normal_Eval;

Global_Search *GlobalSearchMethod = NULL;
Local_Search *LocalSearchMethod = NULL;
//// Local_Search *UnboundLocalSearchMethod = NULL;

//Declaration of Variables for particle swarm optimization (PSO) 
// (defaults now set in constructor, see pso.h)
//PSO in SODOCK 
// MP these max values are as received from ??? and seem excessive...
float pso_tvmax = 2.0;
float pso_qvmax = 1.0;
float pso_rvmax = DegreesToRadians(50.0);
PSO_Options pso_options ; //  MP in progress


info = (GridMapSetInfo *) calloc(1, sizeof(GridMapSetInfo) );
if(info == NULL) stop("failed to allocate grid info structure");
ad_energy_tables = (EnergyTables *) calloc(1, sizeof(EnergyTables) );
if(ad_energy_tables == NULL)  stop("failed to allocate energy tables");
unbound_energy_tables = (EnergyTables *) calloc(1, sizeof(EnergyTables) );
if(unbound_energy_tables == NULL)  stop("failed to allocate unbound energy tables");

// Create a coordinate at the origin:
Coord origin;
origin.x = 0.;
origin.y = 0.;
origin.z = 0.;

//______________________________________________________________________________
/*
** Get the time at the start of the run...
*/

jobStart = times( &tms_jobStart );

#ifdef _OPENMP
/*
** OpenMP initialization
*/

 /* make sure we use no more than compile-time max number of threads,
  * controlled by size of arrays in the per-thread random number generators, com.cc.
  * For example, if NUMG is 8, multi-core or multi-CPU processors up to 8
  * hardware threads are supported.
  */
   if(omp_get_max_threads()>NUMG) omp_set_num_threads(NUMG);
#endif

//_____________________________________________________________________________
/*
** Boinc initialization
*/
#ifdef BOINC
    int flags = 0;
    int rc;
    flags =
      BOINC_DIAG_DUMPCALLSTACKENABLED |
      BOINC_DIAG_HEAPCHECKENABLED |
      BOINC_DIAG_REDIRECTSTDERR |
      BOINC_DIAG_REDIRECTSTDOUT ;
    boinc_init_diagnostics(flags);

#ifdef BOINCCOMPOUND
    BOINC_OPTIONS options;
    options.main_program = false;
    options.check_heartbeat = false; // monitor does check heartbeat
    options.handle_trickle_ups = false;
    options.handle_trickle_downs = false;
    options.handle_process_control = false;
    options.send_status_msgs = true;// only the worker programs (i.e. model) sends status msgs
    options.direct_process_action = true;// monitor handles suspend/quit, but app/model doesn't
    // Initialisation of Boinc
    rc =  boinc_init_options(options); //return 0 for success
    if( rc ){
      fprintf(stderr,"BOINC_ERROR: boinc_init_options() failed \n");
      exit(rc);
    }

#else
    // All BOINC applications must initialize the BOINC interface:
    rc = boinc_init();
    if (rc){
      fprintf(stderr, "BOINC_ERROR: boinc_init() failed.\n");
      exit(rc);
    }
#endif
#endif

// set initial outlev value
(void) getoutlev("default", &outlev); // see bottom of main.cc and constants.h
if(outlev!=LOGFORADT) stop("default outlev fail"); // debug

//______________________________________________________________________________
/*
** Parse the arguments in the command line...
** setflags() conditionally sets globals dock_param_fn, parFile, logFile
** and others
*/

if ( setflags(argc,argv,version_num.c_str()) == -1) {
    exit(EXIT_FAILURE);
} /* END PROGRAM */

// do not allow parameter file reading from standard input (stdin)
// when the standard input is a terminal - almost certainly an error
if ( 0 == strlen(dock_param_fn) && isatty( fileno(stdin)) ) {
	stop("no parameter file (.dpf) specified and AutoDock input is a terminal");
}


//______________________________________________________________________________
/*
** Initialize torsion arrays and constants.
** This should be done for each new ligand read 
*/


int ltorfmt = 4;
(void) strcpy( torfmt, "%*s" ); /* len(torfmt) is 3 chars */

for (j = 0;  j < MAX_ATOMS;  j++ ) {
    type[j] = 0;
    ignore_inter[j] = 0;
}

for (i = 0; i  < MAX_TORS;  i++ ) {
    for (j = 0;  j < MAX_ATOMS;  j++ ) {
        tlist[i][j] = 0;
    }
}

for (i = 0; i  < MAX_TORS;  i++ ) {
    B_isTorConstrained[i] = 0;
    US_torProfile[i][0] = 0;
    N_con[i] = 0;
}

initialiseState( &sInit );
initialiseState( &(ligand.S) );

initialiseQuat( &q_reorient ); // set to identity

F_W      =  360.0 / NTORDIVS;
F_hW     =  F_W  / 2.0;
F_A_from = -360.0 + F_hW;
F_A_to   =  360.0 + F_hW;

for (k = 0; k < MAX_RUNS; k++) {
    for (i = 0; i  < MAX_TORS;  i++ ) {
        sHist[k].tor[i] = 0.0;
    }
}

for (i = 0; i < MAX_TORS;  i++ ) {
    if ( (ltorfmt += 4) > LINE_LEN ) {
        prStr( error_message, "%s:  ERROR: MAX_TORS = %d torsions declared in \"constants.h\";\n\t LINE_LEN = %d, Therefore you must change \"LINE_LEN\" to exceed %d...\n", programname, MAX_TORS, LINE_LEN, 4+4*MAX_TORS );
        stop( error_message );
    } else {
        (void) strcat( torfmt, " %lf" );  /* add on 4 chars  for each new torsion... */
    }
} /* len(torfmt) is 3+4*MAX_TORS chars */

for (j = 0; j < MAX_NONBONDS; j++) {
    nonbondlist[j].a1 = nonbondlist[j].a2 = 0;
}

for (j = 0; j < MAX_RUNS; j++) {
    econf[j] = 0.0;
}

B_constrain_dist = B_haveCharges = FALSE;
ntor = atomC1 = atomC2 = 0;
sqlower = squpper = 0.0;

if (clktck == 0) {        /* fetch clock ticks per second first time */
    if ( (clktck = sysconf(_SC_CLK_TCK)) < (long)0L) {
        stop("\"sysconf(_SC_CLK_TCK)\" command failed in \"main.c\"\n");
    } else {
        idct = (Real)1.0 / (Real)clktck;
        if (debug) {
            pr(logFile, "N.B. debug is on and set to %d\n\n", debug);
            pr(logFile, "\n\nFYI:  Number of clock ticks per second = %d\n", (int)clktck);
            pr(logFile, "FYI:  Elapsed time per clock tick = %.3e milli-seconds\n\n\n\n", idct * 1000. );
        }
    }
}

(void) strcpy(FN_rms_ref_crds,"unspecified filename\0");


//______________________________________________________________________________
/*
** log(x): compute the natural (base e) logarithm of x,
*/

F_lnH = ((Real)log(0.5));



//______________________________________________________________________________
/*
** Determine initial output level before we output anything.
** We must parse the entire DPF -- silently -- for any outlev settings
** or flexible residues file specification
*/

while( fgets(line, LINE_LEN, parFile) != NULL ) { /* Pass 1 PARSING-DPF parFile */
    dpf_keyword = parse_dpf_line( line );
    if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';  // remove newline if last char in line

    switch( dpf_keyword ) {
    case DPF_OUTLEV:
        /*
        **  outlev
        **  Output level,
	**  syntax errors found in this first pass could lead to crypic error
	**  messages since the DPF lines are not echoed in the first pass
	**  
        */
	if(! getoutlev(line, &outlev)) {
	char msg[LINE_LEN+60];
	sprintf(msg, "syntax error or illegal value in DPF outlev setting '%s'", line);
	stop(msg);
	}
        break;

    case DPF_FLEXRES:
        // The DPF specifies a flexible residues file
        // -- set a flag
        // -- get the filename
        B_have_flexible_residues = TRUE;
        get1arg( line, "%*s %s", FN_flexres, "FLEXIBLE_RESIDUES" );
        break;

    default:
        break;
    } // switch( dpf_keyword )
} // while

// Rewind DPF, so we can resume normal parsing
(void) rewind( parFile );


//______________________________________________________________________________
/*
** Output banner, date/time of run, hostname, working directory
*/

banner( version_num.c_str(), outlev, logFile);

if ( outlev >= LOGBASIC ) {
(void) fprintf(logFile, "                     main.cc  $Revision: 1.213 $\n\n");
(void) fprintf(logFile, "                   Compiled on %s at %s\n\n\n", __DATE__, __TIME__);
}

(void) strcpy(hostnm, "unknown_host");
#ifdef HAVE_GETHOSTNAME
gethostname( hostnm, sizeof hostnm );
#endif
if(hostnm[0]=='\0') strcpy(hostnm, "unknown_host");

(void) strcpy(FN_current_working_directory, "unknown_directory");
if(NULL==getcwd(FN_current_working_directory, sizeof FN_current_working_directory))
  strcpy(FN_current_working_directory, "unknown_directory");

if(outlev>=LOGMIN) {
	pr( logFile, "This file was created at:\t\t\t" );
	printdate( logFile, 1 );
	pr( logFile, "                   on host:\t\t\"%s\"\n", hostnm );
	pr(logFile, "Current Working Directory = \"%s\"\n", FN_current_working_directory);
}


//______________________________________________________________________________
if(outlev>=LOGFORADT) {

(void) fprintf(logFile, "\n      ________________________________________________________________\n\n");
(void) fprintf(logFile, "                   SETTING UP DEFAULT PARAMETER LIBRARY\n");
(void) fprintf(logFile, "      ________________________________________________________________\n\n\n");
}

//______________________________________________________________________________
//
// Read in default parameters
//
setup_parameter_library(logFile, outlev, "default Unbound_Same_As_Bound", Unbound_Same_As_Bound, &AD4);

//
// Compute the look-up table for the distance-dependent dielectric function
//
if(outlev >= LOGETABLES)
(void) fprintf(logFile, "\n\nPreparing Energy Tables for Bound Calculation:\n\n");
setup_distdepdiel(logFile, outlev, ad_energy_tables);
if(outlev >= LOGETABLES)
(void) fprintf(logFile, "Preparing Energy Tables for Unbound Calculation:\n\n");
setup_distdepdiel(logFile, outlev, unbound_energy_tables);


// set initial default seeds for random number generator (function is below, at end of main.cc)
set_seeds( seed, seedIsSet, runseed, outlev, logFile);

//______________________________________________________________________________

if(outlev>LOGFORADT) {
(void) fprintf(logFile, "\n      ___________________________________________________\n\n");
(void) fprintf(logFile,   "             PARSING INPUT DOCKING PARAMETER FILE\n");
(void) fprintf(logFile,   "      ___________________________________________________\n\n");
}

//______________________________________________________________________________
/*
** (Note: "dock_param_fn" set in "setflags.c"...)
*/
pr( logFile, "Docking parameter file (DPF) used for this docking:\t\t%s\n", dock_param_fn );

//______________________________________________________________________________
/*
** Start reading in the DPF parameter/run-control file,
*/

while( fgets(line, LINE_LEN, parFile) != NULL ) { /* Pass 2 PARSING-DPF parFile */
    // "line" is a string containing the current line of the input DPF.

    (void) fflush(logFile);
    dpf_keyword = parse_dpf_line( line );

    switch( dpf_keyword ) {
        case -1:
            sprintf( error_message,
               "DPF> %s\n%s: ERROR: Unrecognized keyword in docking parameter file.\n",
               line, programname );
            stop( error_message );

            break;

        case DPF_BLANK_LINE:
        case DPF_COMMENT:
	    if(outlev>=LOGBASIC) pr( logFile, "DPF> %s\n", line );
            break;

        default:
	    if(outlev>LOGFORADT) pr(logFile, "\n\n");
            if(outlev>=LOGBASIC) pr( logFile, "DPF> %s\n", line );
            indcom = strindex( line, "#" );
	    /* truncate line at comment - # mark must be first character
	     * in line, or else preceded by white space (blank or tab)
	     */
            if (indcom != -1 && 
	      (indcom==0 || (isascii(line[indcom-1]) && isspace(line[indcom-1]))) 
	      ) line[ indcom ] = '\0';
            break;
    } /* switch */

    switch( dpf_keyword ) {

//______________________________________________________________________________

    case DPF_BLANK_LINE:
    case DPF_COMMENT:
        break;

//______________________________________________________________________________

    case DPF_PARAMETER_VERSION:
        /*
        ** autodock_parameter_version string
        **
        **
        **
        ** initial implementation ignores value of string
        */

        B_found_autodock_parameter_version = 1==sscanf( line, "%*s %s", autodock_parameter_version );
        pr( logFile, "\tAutodock parameter version %s.\n", autodock_parameter_version );
        break;

/*____________________________________________________________________________*/

    case DPF_OUTPUT_POP_STATS:
        /*
        **  output_population_statistics [option string]
	*  must be after outlev line to have effect
        */
        nfields = sscanf( line, "%*s %s %d %d", 
	  c_mode_str, 
	  &output_pop_stats.everyNgens,
	  &output_pop_stats.everyNevals );
        if (nfields==3 && streq(c_mode_str, "basic")) {
	    output_pop_stats.level = 1; // basic
	    }
	    // nothing besides "basic (int) (int)" is supported yet
        else stop("unsupported option in \"output_population_statistics\" line.\n");
	break;
/*____________________________________________________________________________*/

    case DPF_OUTPUT_RESNUM_AS:
        /*
        ** pdbqt format for residues in dlgs
        *  default is to keep the residue number string from input
        *  possible values: 'resnum' and 'runnum'
        *  intended to replace '-k' option in setflags.cc
        */
        nfields = sscanf( line, "%*s %s", c_mode_str );
        if (nfields==1 && streq(c_mode_str, "resnum")) {
	    keepresnum = TRUE; // default
	    }
        else if (nfields==1 && streq(c_mode_str, "runnum")) {
	    keepresnum = FALSE; // old -k option
	    }
        else stop("unsupported option in \"output_resnum_as\" line.\n");
	break;
/*____________________________________________________________________________*/

    case DPF_OUTLEV:
        /*
        **  outlev
        **  Output level,
        */
	if(! getoutlev(line, &outlev)) syntaxstop("outlev");
        output_pop_stats.everyNgens = (unsigned int) OUTLEV0_GENS; // default
        pr( logFile, "Output Level = %d " , outlev);
        switch ( outlev ) {
	case LOGMIN:
            pr( logFile, "ONLY STATE VARIABLES OUTPUT, NO COORDINATES.\n" );
            output_pop_stats.everyNgens = (unsigned int) OUTLEV0_GENS;
            break;
	case LOGMINCLUST:
            pr( logFile, "ONLY STATE VARIABLES AND CLUSTERING OUTPUT, NO COORDINATES.\n" );
            output_pop_stats.everyNgens = (unsigned int) OUTLEV0_GENS;
            break;
        case LOGBASIC:
            pr( logFile, " BASIC OUTPUT DURING DOCKING (LOGBASIC).\n" );
            output_pop_stats.everyNgens = (unsigned int) OUTLEV0_GENS;
            break;
        case LOGFORADT:
            pr( logFile, " ADT-COMPATIBLE OUTPUT DURING DOCKING.\n" );
            output_pop_stats.everyNgens = (unsigned int) OUTLEV0_GENS;
            break;
        case LOGRUNV:
            pr( logFile, " EXPANDED OUTPUT DURING DOCKING.\n");
            output_pop_stats.everyNgens = (unsigned int) OUTLEV1_GENS;
	    break;
	case LOGLIGREAD:
	    pr(logFile, " EXPANDED OUTPUT DURING LIGAND SETUP.\n" );
	    break;
	case LOGRECREAD:
	    pr(logFile, " EXPANDED OUTPUT DURING LIGAND/RECEPTOR SETUP.\n" );
            break;
        case LOGRUNVV:
        case LOGRUNVVV:
            pr( logFile, " FULL OUTPUT DURING DOCKING.\n");
            output_pop_stats.everyNgens = (unsigned int) OUTLEV2_GENS;
	    break;
	case LOGETABLES:
	case LOGNBINTE:
	case LOGNBINTEV:
	    pr(logFile, " EXPANDED OUTPUT DURING ENERGY TABLE SETUP.\n" );
	    break;
	default:
	    pr(logFile, " WARNING, undocumented outlev setting %d.\n", outlev );
	    break;
        }
        if(output_pop_stats.everyNgens>0) pr( logFile, "\n\tOutput population statistics every %u generations.\n", output_pop_stats.everyNgens );
        else if(outlev>LOGFORADT) pr( logFile, "\n\tNever output generation-based population statistics.\n");
        if(output_pop_stats.everyNevals>0) pr( logFile, "\n\tOutput population statistics every %u energy evaluations.\n", output_pop_stats.everyNevals );
        else if(outlev>LOGFORADT) pr( logFile, "\n\tNever output evaluation-count-based population statistics.\n");
        break;

/*____________________________________________________________________________*/

    case DPF_PARAMETER_LIBRARY:
        /*
        ** parameter_file AD4_parameters.dat
        **  or
        ** parameter_library AD4_parameters.dat
        **
        ** initial implementation based on hsearch was suggested by Mike Pique
        */

        parameter_library_found = 1==sscanf( line, "%*s %s", FN_parameter_library );
        read_parameter_library(logFile, outlev, FN_parameter_library, &AD4);

        break;

/*____________________________________________________________________________*/

    case DPF_INCLUDE_1_4_INTERACTIONS:
        /*
         * include_1_4_interactions 0.5
         *
         * Set the Boolean variable, B_include_1_4_interactions, to TRUE.
         *
         * NOTE:  You must use this command _before_ the "move ligand.pdbqt"
         *        command, since "include_1_4_interactions" affects how the Ligand
         *        PDBQT specified by the "move" command will be interpreted.
         */
        if (B_found_move_keyword == TRUE) {
            // If we have found the move keyword already, warn the user
            // that this command ("include_1_4_interactions 0.5") should have
            // been given before this!
            pr(logFile, "this INCLUDE_1_4_INTERACTIONS command must be before the \"move ligand.pdbqt\" command, since this command affects how the PDBQT file will be interpreted.\n\n");
	    stop("");
        }
        get1arg( line, "%*s " FDFMT, &scale_1_4 , "INCLUDE_1_4_INTERACTIONS");
        B_include_1_4_interactions = TRUE;
        print_1_4_message(B_include_1_4_interactions, scale_1_4, outlev, logFile);
        break;

//______________________________________________________________________________


    case DPF_SCALE_EINTERMOL:
        /*
        **  scale_eintermol
        **  re-scale intermolecular energy term
        */
        get1arg( line, "%*s " FDFMT, &scale_eintermol, "SCALE_EINTERMOL");
	  pr(logFile,"  Intermolecular energy term will be scaled by factor %f\n", scale_eintermol);
	 break;
//______________________________________________________________________________


    case DPF_SMOOTH:
        /*
        **  smooth r_smooth
        **  set internal non-bond table smoothing range (not radius) in Angstroms
        ** Typical value of r_smooth is 0.5 Angstroms
        */
	if(B_found_ligand_types) {
		stop("You must specify the smoothing range before specifying the ligand_types.");
		}
        get1arg( line, "%*s " FDFMT, &r_smooth, "SMOOTH");
        (void) fprintf( logFile, "\nInternal energy non-bond potentials will be smoothed over range %.3lf Angstrom\n\n", r_smooth);
	 break;
//______________________________________________________________________________


    case DPF_INTELEC:
        /*
        **  intelec  [ off ]
        **  Calculate internal electrostatic energies...
        */
        nfields = sscanf( line, "%*s %s", param[0]);
	if ( streq(param[0], "off")) {
	   B_calcIntElec = FALSE;
	   if (outlev >= LOGBASIC) pr( logFile,
	   "Electrostatic energies will not be calculated for interactions between moving atoms.\n");
	} else {
	   B_calcIntElec = TRUE;
	   if (outlev >= LOGBASIC) pr( logFile, 
	   "Electrostatic energies will be calculated for all non-bonds between moving atoms.\n");
        }
	// check for obsolete numeric parameter
        nfields = sscanf( line, "%*s " FDFMT, &AD3_FE_coeff_estat );
        if (nfields == 1) {
                pr(logFile, "NOTE!  Internal electrostatics will NOT be scaled by the factor specified by this command,  %.4f -- the coefficient set by this command is ignored in AutoDock 4;\n", AD3_FE_coeff_estat);
                pr(logFile, "       the coefficient that will actually be used should be set in the parameter library file.\n");
                pr(logFile, "       The coefficient for the electrostatic energy term is %.4f", AD4.coeff_estat);
                if (parameter_library_found) {
                    pr( logFile, " as specified in parameter library \"%s\".\n", FN_parameter_library );
                } else {
                    pr( logFile, ", the factory default value.\n");
                }

	    stop("illegal obsolete numeric value in intelec line");
            }
        break;

//______________________________________________________________________________

    case DPF_SEED:
        /*
        **  seed
        **  Set the random-number generator's seed value,
        */
        nfields = sscanf( line, "%*s %s %s", param[0], param[1]);
        //pr(logFile, "%d seed%s found.\n", nfields, pl(nfields));
        if ((nfields==2) || (nfields==1)) {
            for (i=0; i<nfields ; i++ ) {
                if (streq(param[i], "time")||streq(param[i],"tim")) {
		    time_t time_seed;
                    seedIsSet[i] = 'T';
		    do {
		        /* seeds<=1 are invalid */
                        seed[i] = (FourByteLong)time( &time_seed );
		    } while ( seed[i]<=1 );
		
		    if(outlev>=LOGRUNV)
                    pr(logFile,"Random number generator seed %d was seeded with the current time, value = " FBL_FMT "\n",i,seed[i]);
                } else if (streq(param[i], "pid")) {
                    seedIsSet[i] = 'P';
                    seed[i] = processid();
		    if(outlev>=LOGRUNV)
                    pr(logFile,"Random number generator seed %d was seeded with the process ID, value   = " FBL_FMT "\n",i,seed[i]);
                } else {
                    seedIsSet[i] = 'U';
                    seed[i] = atol(param[i]);
		    if(seed[i]<=1) stop("Random number seed cannot be zero or one, or negative");
		    if(outlev>=LOGRUNV)
                    pr(logFile,"Random number generator seed %d was seeded with the user-specified value  " FBL_FMT "\n",i,seed[i]);
                }
            }/*i*/
	set_seeds( seed, seedIsSet, runseed, outlev, logFile);
        } else stop("Error encountered reading SEED line");

	
	/* debugging extension: if a third field is present, write out that
	 * many random numbers (as integers) to a private file then exit AutoDock
	 */
        nfields = sscanf( line, "%*s %*s %*s %d", &i);
	if(nfields==1) {
		FILE *rand_fd;
		if(NULL!=(rand_fd=fopen("randoms", "w")) ) {
			int j;
			for(j=0;j<i;j++) fprintf(rand_fd, FBL_FMT "\n", ignlgi());
			fclose(rand_fd);
			exit(0);
			}
		}

        break;

/*____________________________________________________________________________*/

    case DPF_LIGAND_TYPES:
        /*
         *  Read in the ligand atom type names, e.g.
         *
         *  ligand_types C HD OA P               # ligand atom type names
         *
         *  The order of the arguments is the index that will
         *  be used for look up in the grid maps, "map_index".
         */

        //  Use the function "parsetypes" to read in the atom types;
        //
        //  The array "ligand_atom_type_ptrs" is returned, having been filled with pointers
        //  to the beginning of each "atom type word" (not atom type characters);
        //  In AutoDock 4, an atom type can be either 1 or 2 characters long.
	//  caution: "parsetypes" modifies its first (string) argument
        //
        num_atom_types = parsetypes(line, ligand_atom_type_ptrs, MAX_ATOM_TYPES);
        if (num_atom_types<0){
            prStr( error_message, "%s:  ERROR! Too many atom types have been found: maximum is %d; we cannot continue !\n\n", programname, MAX_ATOM_TYPES );
            pr_2x( logFile, stderr, error_message );
	    stop(error_message);
        }



        B_found_ligand_types = TRUE;
        info->num_atom_types = num_atom_types;

        for (i=0; i<num_atom_types; i++) {
            strcpy(info->atom_type_name[i], ligand_atom_type_ptrs[i]);
            if (!strncmp(&info->atom_type_name[i][0], "H", 1)) h_index = i;
#ifdef DEBUG
            (void) fprintf(logFile, "%d %s ->%s\n",i, ligand_atom_type_ptrs[i], info->atom_type_name[i]);
            (void) fprintf(logFile, "h_index =%d\n",h_index);
#endif
        }

        if (num_atom_types > 0) {
            B_atom_types_found = TRUE;
        } else {
            prStr( error_message, "%s:  ERROR! No atom types have been found; we cannot continue without this information!\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            prStr( error_message, "%s:  ERROR! Are you trying to use an AutoDock 3 DPF with AutoDock 4?\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            stop(error_message);
        }

        if (debug > 0) {
            for (i=0; i<num_atom_types; i++) {
                (void) fprintf(logFile, "info->atom_type_name[%d] = %s\n", i, info->atom_type_name[i] );
            }
        }

        // For all ligand atom types... set up the map_index
        // "ligand_types"
        for (i=0; i<num_atom_types; i++) {
            foundParameter = apm_find(info->atom_type_name[i]);
            if (foundParameter != NULL ) {
                // Not NULL means we have found this atom type's parameters.
                // Set the ParameterEntry's "map_index" member to the
                // 0-based index it had in the list of ligand types supplied in the DPF "types" line:
                foundParameter->map_index = i;
                parameterArray[i] = *(foundParameter);
                if (outlev >= LOGLIGREAD) {
                    (void) fprintf( logFile, 
		    "Parameters found for ligand type \"%s\" (grid map index = %d, weighted well depth, epsilon = %6.4f, Rij = %6.4f)", 
		    foundParameter->autogrid_type, foundParameter->map_index, 
		    foundParameter->epsij, foundParameter->Rij);
                    if (parameter_library_found) {
                        pr( logFile, " in parameter library \"%s\".\n", FN_parameter_library );
                    } else {
                        pr( logFile, "\n");
                    }
                }
            } else {
                // We could not find this parameter -- return error here
                prStr( error_message,"%s: ERROR:  Unknown ligand atom type \"%s\"; add parameters for it to the parameter library first!\n", programname, info->atom_type_name[i]);
                pr_2x( logFile, stderr, error_message );
                if (parameter_library_found) {
                    prStr( error_message,"%s:         Edit the parameter library file \"%s\" and try again.\n", programname, FN_parameter_library );
                    pr_2x( logFile, stderr, error_message );
                }
                stop(error_message);
            } // if / else apm_find
        } // for i

        // Calculate the internal energy table

        // loop over atom types, i
        for (i=0; i<num_atom_types; i++) {

            //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
            //  Lennard-Jones and Hydrogen Bond Potentials

	    // i
	    double Ri, epsi, Ri_hb, epsi_hb;
	    hbond_type hbondi;

            Ri = parameterArray[i].Rij;
            epsi = parameterArray[i].epsij;
            Ri_hb = parameterArray[i].Rij_hb;
            epsi_hb = parameterArray[i].epsij_hb;
            hbondi = parameterArray[i].hbond;

            // loop over atom types, j, from i to number of atom types
            for (j=i; j<num_atom_types; j++) {

		// j
		double Rj, epsj, Rj_hb, epsj_hb, epsij, Rij;
		hbond_type hbondj;
		Boole is_hbond;
                //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
                //  Lennard-Jones and Hydrogen Bond Potentials

                Rj = parameterArray[j].Rij;
                epsj = parameterArray[j].epsij;
                Rj_hb = parameterArray[j].Rij_hb;
                epsj_hb = parameterArray[j].epsij_hb;
                hbondj = parameterArray[j].hbond;

                // we need to determine the correct xA and xB exponents
                xA = 12; // for both LJ, 12-6 and HB, 12-10, xA is 12
                xB =  6; // assume we have LJ, 12-6

                if ( ((hbondi == DS) || (hbondi == D1)) && ((hbondj == AS) || (hbondj == A1) || (hbondj == A2)) ) {
                    // i is a donor and j is an acceptor.
                    // i is a hydrogen, j is a heteroatom
                    Rij = Rj_hb;
                    epsij = epsj_hb;
                    xB = 10;
		    is_hbond = TRUE;
                } else if ( ((hbondi == AS) || (hbondi == A1) || (hbondi == A2)) && ((hbondj == DS) || (hbondj == D1))) {
                    // i is an acceptor and j is a donor.
                    // i is a heteroatom, j is a hydrogen
                    Rij = Ri_hb;
                    epsij = epsi_hb;
                    xB = 10;
		    is_hbond = TRUE;
                } else {
                    // we need to calculate the arithmetic mean of Ri and Rj
                    Rij = arithmetic_mean(Ri, Rj);
                    // we need to calculate the geometric mean of epsi and epsj
                    epsij = geometric_mean(epsi, epsj);
		    is_hbond = FALSE;
                }

                /* Check that the Rij is reasonable */
                if ((Rij < RIJ_MIN) || (Rij > RIJ_MAX)) {
                    (void) fprintf( logFile,
                    "WARNING: pairwise distance, Rij, %.2f, is not a very reasonable value for the equilibrium separation of two atoms! (%.2f Angstroms <= Rij <= %.2f Angstroms)\n\n", Rij, RIJ_MIN, RIJ_MAX);
                    (void) fprintf( logFile, "Perhaps you meant to use \"intnbp_coeffs\" instead of \"intnbp_r_eps\"?\n\n");
                    /* gmm commented out for dave goodsell, mutable atoms
                     * exit(EXIT_FAILURE); */
                }
                /* Check that the epsij is reasonable */
                if ((epsij < EPSIJ_MIN) || (epsij > EPSIJ_MAX)) {
                    (void) fprintf( logFile,
                    "WARNING: well-depth, epsilon_ij, %.2f, is not a very reasonable value for the equilibrium potential energy of two atoms! (%.2f kcal/mol <= epsilon_ij <= %.2f kcal/mol)\n\n", epsij, EPSIJ_MIN, EPSIJ_MAX);
                    (void) fprintf( logFile, "Perhaps you meant to use \"intnbp_coeffs\" instead of \"intnbp_r_eps\"?\n\n");
                    /* gmm commented out for dave goodsell, mutable atoms
                     * exit(EXIT_FAILURE); */
                }
                /* Defend against division by zero... */
                if (xA != xB) {
		    double tmpconst = epsij / (Real)(xA - xB);
                    cA = tmpconst * pow( (double)Rij, (double)xA ) * (Real)xB;
                    cB = tmpconst * pow( (double)Rij, (double)xB ) * (Real)xA;

		    if(outlev >= LOGETABLES) {
                    pr(logFile, "\nCalculating internal non-bonded interaction energies for docking calculation;");
                    pr(logFile, "\n smoothing range is %.4f (i.e., d - %2f to d + %.2f Angstrom)\n", 
		     r_smooth, r_smooth/2, r_smooth/2);
		    }
                    intnbtable( &B_havenbp, a1, a2, info, cA, cB, xA, xB, is_hbond, 
		      r_smooth, AD4, sigma, ad_energy_tables, BOUND_CALCULATION, 
		      logFile, outlev);
		    if(outlev>=LOGETABLES)
                    pr(logFile, "\nCalculating internal non-bonded interaction energies for unbound conformation calculation;\n");
                    intnbtable( &B_havenbp, a1, a2, info, cA_unbound, cB_unbound, xA_unbound, xB_unbound, is_hbond,
		      r_smooth, AD4,  sigma, unbound_energy_tables, UNBOUND_CALCULATION, 
		      logFile, outlev);
                    // Increment the atom type numbers, a1 and a2, for the internal non-bond table
                    a2++;
                    if (a2 >= info->num_atom_types) {
                        a1++;
                        a2 = a1;
                    }

                } else {
                    pr(logFile,"ERROR: Exponents must be different, to avoid division by zero!\n\tAborting...\n");
		    stop("exponent would cause division by zero");
                }

            } // for j
        } // for i
        break;

//______________________________________________________________________________

    case DPF_FLD:
        /*
        ** fld
        ** GRID_DATA_FILE
        ** Read the (AVS-format) grid data file, .fld
	**
	** Fatal error if "ligand_types" has not already appeared.
        */
	if(! B_found_ligand_types) {
		stop("You must specify the ligand_types before reading the grid data file.");
		}
        readfield( info, line, jobStart, tms_jobStart, outlev, logFile );
        num_maps = 0;

        // Dynamically allocate memory for the maps, clear to 0.
	// We need space for all the atom maps (info->num_atom_types), 
	// plus the electrostatic potential and the desolvation map

	free(map); // note: it is OK to free even if NULL
	info->num_all_maps = info->num_atom_types+NUM_NON_VDW_MAPS;
	info->num_alloc_maps = info->num_all_maps;
	info->map_alloc_size =
	   info->num_alloc[Z] * info->num_alloc[Y] * 
	   info->num_alloc[X] * info->num_alloc_maps;
        map = (MapType *) calloc(info->map_alloc_size, sizeof(MapType));
	if(outlev >= LOGRECREAD) {
	   pr(logFile, "Allocating %d x %d x %d (x,y,z) grid of %d maps, %ld bytes\n", 
	   info->num_alloc[X] , info->num_alloc[Y] , 
	   info->num_alloc[Z] , info->num_alloc_maps,
	   (long)info->map_alloc_size * sizeof(MapType));
	   }

        if (map == NULL) {
	   pr(logFile, "Failed to allocate %d x %d x %d (x,y,z) grid of %d maps, %ld bytes\n", 
	   info->num_alloc[X] , info->num_alloc[Y] , 
	   info->num_alloc[Z] , info->num_alloc_maps,
	   (long)info->map_alloc_size * sizeof(MapType));
            prStr(error_message, "%s:  Sorry, there is not enough memory to store the grid maps.  Please use smaller maps and/or fewer atom types.\n", programname);
            stop(error_message);
        }
        break;

//______________________________________________________________________________

    case DPF_MAP:
        /*
        ** map
        ** ATOMIC AFFINITY, ELECTROSTATIC POTENTIAL OR DESOLVATION ENERGY GRID MAP
        ** Read in active site grid map...
        */
        B_charMap = FALSE;
        if (B_atom_types_found == TRUE) {
            // Read in the AutoGrid atomic affinity map
            // map_index could be incremented here if we had the atom_type stored in each map...
            map_stats = readmap( line, outlev, jobStart, tms_jobStart, B_charMap, &B_havemap, num_maps, info, map, 'a', logFile);
            if( outlev >= LOGRECREAD ) pr(logFile, "Min= %.3lf Mean= %.3lf Max= %.3lf\n\n",
                    map_stats.minimum, map_stats.mean, map_stats.maximum);
            num_maps++;

        } else {
            prStr( error_message, "%s:  ERROR! No atom types have been found; we cannot continue without this information!\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            prStr( error_message, "%s:  ERROR! Are you trying to use an AutoDock 3 DPF with AutoDock 4?\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            stop(error_message);
        }
        break;

//______________________________________________________________________________

    case DPF_ELECMAP:
        /*
         *  elecmap file.e.map
         */
        map_stats = readmap( line, outlev, jobStart, tms_jobStart, B_charMap, &B_havemap, num_maps, info, map, 'e', logFile);
        if( outlev >= LOGRECREAD ) pr(logFile, "Min= %.3lf Mean= %.3lf Max= %.3lf\n\n",
                map_stats.minimum, map_stats.mean, map_stats.maximum);
        ElecMap = num_maps;
        B_found_elecmap = TRUE;
        num_maps++;
        break;

//______________________________________________________________________________

    case DPF_DESOLVMAP:
        /*
         *  desolvmap file.d.map
         */
        map_stats = readmap( line, outlev, jobStart, tms_jobStart, B_charMap, &B_havemap, num_maps, info, map, 'd', logFile);
        if( outlev >= LOGRECREAD ) pr(logFile, "Min= %.3lf Mean= %.3lf Max= %.3lf\n\n",
                map_stats.minimum, map_stats.mean, map_stats.maximum);
        DesolvMap = num_maps;
        B_found_desolvmap = TRUE;
        num_maps++;
        break;

//______________________________________________________________________________

    case DPF_CHARMAP:
        /*
        ** charmap
        ** ATOMIC AFFINITY, ELECTROSTATIC POTENTIAL OR DESOLVATION ENERGY GRID MAP
        ** Read in active site grid map...
        */
        B_charMap = TRUE;
        if (B_atom_types_found == TRUE) {
            // map_index could be incremented here if we had the atom_type stored in each map...
            map_stats = readmap( line, outlev, jobStart, tms_jobStart, B_charMap, &B_havemap, num_maps, info, map, 'c', logFile);
            if( outlev >= LOGRECREAD ) pr(logFile, "Min= %.3lf Mean= %.3lf Max= %.3lf\n\n",
                    map_stats.minimum, map_stats.mean, map_stats.maximum);
            num_maps++;
        } else {
            prStr( error_message, "%s:  ERROR! No atom types have been found; we cannot continue without this information!\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            prStr( error_message, "%s:  ERROR! Are you trying to use an AutoDock 3 DPF with AutoDock 4?\n\n", programname );
            pr_2x( logFile, stderr, error_message );
	    stop(error_message);
        }
        break;

//______________________________________________________________________________

    case DPF_MOVE:
        /*
        ** move ligand_file.pdbqt
        ** Specify the movable ligand,
        */
        //
        // Initialisations that must be done before reading in a new ligand...
        //
        if (num_maps != num_atom_types + NUM_NON_VDW_MAPS) { //  dsolv map and elec map
            prStr(error_message, "\n\nMISSING MAP ERROR:\nnumber of maps %d does not match number expected for %d ligand types. \nUnable to continue.\n", num_maps, num_atom_types);
            stop(error_message);
        }

        nconf = 0;
        for (k = 0; k < MAX_RUNS; k++) {
            for (i = 0; i  < MAX_TORS;  i++ ) {
                sHist[k].tor[i] = 0.0;
            }
            econf[k] = 0.0;
        }
        for (j = 0;  j < MAX_ATOMS;  j++ ) {
            type[j] = 0;
            ignore_inter[j] = 0;
        }
        for (i = 0; i  < MAX_TORS;  i++ ) {
            for (j = 0;  j < MAX_ATOMS;  j++ ) {
                tlist[i][j] = 0;
            }
            B_isTorConstrained[i] = 0;
            US_torProfile[i][0] = 0;
            N_con[i] = 0;
        }
        for (j = 0; j < MAX_NONBONDS; j++) {
            nonbondlist[j].a1 = nonbondlist[j].a2 = 0;
        }
        for (j=0; j<3; j++) {
            Nnb_array[j] = 0;
        }
        initialiseState( &sInit );
        initialiseState( &(ligand.S) );
        initialiseQuat( &q_reorient ); // set to identity
        B_constrain_dist = B_haveCharges = FALSE;
        ntor = atomC1 = atomC2 = 0;
        ntor_ligand = 0;
        ntorsdof = 0;
        sqlower = squpper = 0.0;
        strcpy( FN_pop_file, "");  // means don't print pop_file
        Nnb = 0;
        ligand_is_inhibitor = 1;
        initialise_energy_breakdown(&eb, 0, 0);
        //
        // end of initialization
        //

        // this is the DPF_MOVE section...
        B_found_move_keyword = TRUE;
        B_found_about_keyword = FALSE; //set false by 'move', set true by 'about'
        B_found_tran0_keyword = FALSE;


        print_1_4_message(B_include_1_4_interactions, scale_1_4, outlev, logFile);

        natom=0;
        ligand = readPDBQT( line,
                            num_atom_types,
                            &natom,
                            crdpdb, charge, &B_haveCharges,
                            type, bond_index,
                            pdbaname, FN_ligand, FN_flexres, 
			    B_have_flexible_residues, atomstuff, 
			    &n_heavy_atoms_in_ligand, &true_ligand_atoms,
                            &B_constrain_dist, &atomC1, &atomC2,
                            &sqlower, &squpper,
                            &ntor, &ntor_ligand,
                            tlist, vt,
                            &Nnb, Nnb_array, nonbondlist,
                            jobStart, tms_jobStart, hostnm, &ntorsdof, 
                            ignore_inter,
                            B_include_1_4_interactions,
                            PDBQT_record, end_of_branch,
			    debug, outlev, logFile);

#ifdef DEBUGTLIST 
// MP 2012-04
	for(int t=0;t<=ntor;t++) {
	fprintf(logFile, " tlist[%2d] = %3d %3d %3d : ", t, tlist[t][0]+1, tlist[t][1]+1, tlist[t][2]);
	for(int aii=0;aii<tlist[t][2];aii++) fprintf(logFile,"%2d ",tlist[t][3+aii]+1);
	fprintf(logFile, "\n");
	}
#endif
	// save crdpdb coords as crdorig 
	for(int a=0;a<natom;a++) for(xyz=0;xyz<SPACE;xyz++)  
	  crdorig[a][xyz]=crdpdb[a][xyz];

        // pre-calculate some values we will need later in computing the desolvation energy
        //
        for (i=0; i<natom; i++) {
            abs_charge[i] = fabs(charge[i]);
            qsp_abs_charge[i] = qsolpar * abs_charge[i];
        }
        pr(logFile, "Number of atoms in ligand:  %d\n\n", true_ligand_atoms);
        pr(logFile, "Number of non-hydrogen atoms in ligand:  %d\n\n", n_heavy_atoms_in_ligand);

        pr(logFile, "Number of vibrational degrees of freedom of ligand:  %d\n\n\n", (3 * true_ligand_atoms) - 6 );
        pr(logFile, "Number of torsional degrees of freedom = %d\n", ntorsdof);

        torsFreeEnergy = (Real)ntorsdof * AD4.coeff_tors;

        pr(logFile, "Estimated loss of torsional free energy upon binding = %+.4f kcal/mol\n\n", torsFreeEnergy);

        for (i=0;i<natom;i++) {
            if (ignore_inter[i] == 1 && outlev>=LOGLIGREAD) {
                pr(logFile, "Special Boundary Conditions:\n");
                pr(logFile, "____________________________\n\n");
                pr(logFile, "AutoDock will ignore the following atoms in the input PDBQT file \nin intermolecular energy calculations:\n");
                pr(logFile, "\n(This is because these residue atoms are at the boundary between \nflexible and rigid, and since they cannot move \nthey will not affect the total energy.)\n\n");
                break;
            }
        }
        for (i=0;i<natom;i++) {
            if (ignore_inter[i] == 1 && outlev>=LOGLIGREAD) {
                pr(logFile, "Atom number %d:  %s\n", i+1, atomstuff[i] );
            }
        }
        pr(logFile, "\n");

        if (!B_haveCharges) {
            pr( logFile, "%s: WARNING! No partial atomic charges have been supplied yet.\n\n",programname);
        } else {
            if (Nnb > 0) {
            if (outlev >= LOGLIGREAD) {
               pr(logFile,"Calculating the product of the partial atomic charges, q1*q2, for all %d non-bonded pairs...\n", Nnb);
               pr(logFile," -- Scaled q1*q2 means multiplied by both  %.1lf (for conversion later on to kcal/mol)\n", (double)ELECSCALE);
               pr(logFile,"    and by the AD4 FF electrostatics coefficient, %.4lf\n\n", (double)AD4.coeff_estat);
                pr(logFile,"Non-bonded                           Scaled\n");
                pr(logFile,"   Pair     Atom1-Atom2    q1*q2      q1*q2\n");
                pr(logFile,"__________  ___________  _________  _________\n");
            } //outlev 
            for (i = 0;  i < Nnb;  i++) {
                atm1 = nonbondlist[i].a1;
                atm2 = nonbondlist[i].a2;
                int t1 = nonbondlist[i].t1;
                int t2 = nonbondlist[i].t2;
                nonbondlist[i].desolv =
                       ( parameterArray[t2].vol * (parameterArray[t1].solpar + qsp_abs_charge[atm1])
                       + parameterArray[t1].vol * (parameterArray[t2].solpar + qsp_abs_charge[atm2]) );
                nonbondlist[i].q1q2 = charge[atm1] * charge[atm2];
		nonbondlist[i].is_hbond = ad_energy_tables->is_hbond[t1][t2];  // MPique untested
                if (outlev >= LOGLIGREAD) {
                    pr(logFile,"   %4d     %5d-%-5d  %7.4f",i+1,atm1+1,atm2+1,nonbondlist[i].q1q2);
                }//outlev
                nonbondlist[i].q1q2 *= ELECSCALE * AD4.coeff_estat;
                if (outlev >= LOGLIGREAD) {
                    pr(logFile,"   %8.4f\n",nonbondlist[i].q1q2);
                }//outlev
            } // for
            } // if NNb > 0
        } // else

        sInit.ntor = ligand.S.ntor;
        ++nmol;
        ++nlig;

        break;

/*____________________________________________________________________________*/

    case DPF_FLEXRES:
        /*
         * flexible_residues file.pdbqt
	 * This token is handled in pass 1, above.
         */
        pr(logFile, "The flexible residues will be read in from \"%s\".\n", FN_flexres);
        break;


#ifdef USING_COLINY
/*____________________________________________________________________________*/

    case DPF_COLINY:
    {
        //ostdiostream fstr(logFile);
        //ostdiostream fstr(logFile->_file);
        //CommonIO::set_streams(&fstr,&fstr,&cin);

        struct tms tms_colinyStart;
        struct tms tms_colinyEnd;

        Clock  colinyStart;
        Clock  colinyEnd;

        int coliny_seed;
        char algname[LINE_LEN];
        char nruns_str[LINE_LEN];
        (void) sscanf(line, "%*s %s %d", algname, &nruns);
        (void) sscanf(line, "%*s %s %s", algname, nruns_str);

        if (streq(algname,"help")) {
            std::vector<double> initvec;
            coliny_init(algname, "", 0);
            prStr(error_message, "%s:  ERROR:  no optimizer type specified.", programname);
            stop(error_message);
            exit(EXIT_FAILURE);
        }
        else if (streq(nruns_str,"help")) {
            std::vector<double> initvec;
            coliny_init(algname, nruns_str, 0);
            prStr(error_message, "%s:  ERROR:  no optimizer type specified.", programname);
            stop(error_message);
            exit(EXIT_FAILURE);
        }


            if (nruns+nconf>MAX_RUNS) {
                prStr(error_message, "%s:  ERROR:  %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", 
		programname, nruns+nconf, MAX_RUNS);
                stop(error_message);
                exit(EXIT_FAILURE);
            }
            exit_if_missing_elecmap_desolvmap_about("coliny");

            evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom,
                            info, map, peratomE, nonbondlist, ad_energy_tables,
                            Nnb, Nnb_array, &group_energy,
			    B_calcIntElec, B_isGaussTorCon, B_isTorConstrained, B_ShowTorE,
                            US_TorE, US_torProfile,
                            vt, tlist,
                            crdpdb, sInit, ligand, ignore_inter, B_include_1_4_interactions, scale_1_4, scale_eintermol,
                            unbound_internal_FE,
                            B_use_non_bond_cutoff, B_have_flexible_residues,
			    true_ligand_atoms, outlev, logFile);

            evaluate.compute_intermol_energy(TRUE);

            char domain[1024];
            // NOTE: Coliny enforces the bound constraints, but since the
            // torsion angles are periodic, we simply prevent the optimizer
            // from going too far.
            if (sInit.ntor > 0) {
                sprintf(domain,"[%f,%f] [%f,%f] [%f,%f] [-1000.0,1000.0]^3 [-3.1416,3.1416] [-3.1416,3.1416]^%d",(double)info->lo[X], (double)info->hi[X], (double)info->lo[Y], (double)info->hi[Y], (double)info->lo[Z], (double)info->hi[Z], sInit.ntor);
                //sprintf(domain,"[%f,%f] [%f,%f] [%f,%f] [-1.0,1.1]^3 [-6.2832,12.5664] [-6.2832,12.5664]^%d",(double)info->lo[X], (double)info->hi[X], (double)info->lo[Y], (double)info->hi[Y], (double)info->lo[Z], (double)info->hi[Z], sInit.ntor);
            } else {
                sprintf(domain,"[%f,%f] [%f,%f] [%f,%f] [-1000.0,1000.0]^3 [-3.1416,3.1416]",(double)info->lo[X], (double)info->hi[X], (double)info->lo[Y], (double)info->hi[Y], (double)info->lo[Z], (double)info->hi[Z]);
            }
            pr(logFile, "Number of Coliny %s dockings = %d run%s\n", algname, nruns, pl(nruns))
            pr(logFile, "Search Domain: %s\n", domain);

            //
            // COLINY-SPECIFIC LOGIC - BEGIN
            //

            try {

                std::vector<double> initvec, finalpt;
                // set up initial point
                initvec.resize(7+sInit.ntor);
                initvec[0] = sInit.T.x;
                initvec[1] = sInit.T.y;
                initvec[2] = sInit.T.z;
                initvec[3] = sInit.Q.x;
                initvec[4] = sInit.Q.y;
                initvec[5] = sInit.Q.z;
                initvec[6] = sInit.Q.w;
                for (j=0; j < sInit.ntor ; j++) {
                  initvec[j+7] = DegreesToRadians(sInit.tor[j]);
                }
                coliny_init(algname, domain, sInit.ntor+7);

                for (j=nconf; j<nruns; j++) {
		  Real eintra = 0.0;  // sum of intramolecular energy for the ligand plus that of the protein
		  Real einter = 0.0; // intermolecular energy between the ligand and the protein
                  fprintf( logFile, "\n\tBEGINNING Coliny %s DOCKING\n",algname);
                  pr(logFile, "\nDoing %s run:  %d/%d.\n", algname, j+1, nruns);

                  //coliny uses a single seed
                  coliny_seed = runseed[j][0]+runseed[j][1];
                  pr(logFile, "Seed: %d ["FBL_FMT"+"FBL_FMT"]\n", coliny_seed, runseed[j][0], runseed[j][1], j);
                  (void) fflush(logFile);

                  colinyStart = times(&tms_colinyStart);

                  finalpt.resize( initvec.size() );
                  int neval, niters;
                  coliny_minimize( coliny_seed, initvec, finalpt, neval, niters );
                  //fstr.flush();

				  // get state after this coliny_minimize run 
                  make_state_from_rep( (double *)&(finalpt[0]), int(finalpt.size()), &sHist[nconf], outlev, logFile);

                  pr(logFile, "\nTotal Num Evals: %d\n", neval);
                  sHist[nconf].Center.x = lig_center[X];
                  sHist[nconf].Center.y = lig_center[Y];
                  sHist[nconf].Center.z = lig_center[Z];
                  printState(logFile, sHist[nconf], 2);

                  colinyEnd = times(&tms_colinyEnd);
                  pr(logFile, "Time taken for this %s run:\n", algname);
                  timesyshms(colinyEnd-colinyStart, &tms_colinyStart, &tms_colinyEnd, logFile);
                  pr(logFile, "\n");

                  pr(logFile, "Total number of Energy Evaluations: %d\n", (int)evaluate.evals() );
                  //pr(logFile, "Total number of Iterations:        %d\n", (int)niters);

                  pr(logFile, "\nFinal docked state:\n");
                  pr( logFile, UnderLine );
                  pr( logFile, "\n\n\tFINAL Coliny %s DOCKED STATE\n",algname );
                  pr( logFile,     "\t____________________________________\n\n\n" );
                  (void) fflush(logFile);

                  writePDBQT( j, runseed[nconf], FN_ligand, dock_param_fn, lig_center,
                              sHist[nconf], ntor, &eintra, &einter, natom, atomstuff,
                              crd, peratomE,
                              charge, abs_charge, qsp_abs_charge,
                              ligand_is_inhibitor,
                              torsFreeEnergy,
                              vt, tlist, crdpdb, nonbondlist,
                              ad_energy_tables,
                              type, 
			      Nnb, Nnb_array, &group_energy, true_ligand_atoms,
			      B_calcIntElec,
                              map,
                              ignore_inter,
                              B_include_1_4_interactions, scale_1_4, parameterArray, unbound_internal_FE,
                              info, DOCKED, PDBQT_record, 
			      B_use_non_bond_cutoff, B_have_flexible_residues, ad4_unbound_model,
			      outlev, logFile);

                  // See also "calculateEnergies.cc", switch(ad4_unbound_model)
                  if (ad4_unbound_model == Unbound_Same_As_Bound) {
                      // Update the unbound internal energy, setting it to the current internal energy
                      unbound_internal_FE = eintra;
                  }
                  econf[nconf] = eintra + einter + torsFreeEnergy - unbound_internal_FE;
                  evaluate.reset();

                  ++nconf;

                } // Next run
                if(write_stateFile){
                  fprintf(stateFile,"\t</runs>\n");
                  (void) fflush(stateFile);
                }
                (void) fflush(logFile);
            }
            catch (std::exception& err) {
              (void)fprintf(logFile, "Caught Exception: %s\n", err.what());
              exit(EXIT_FAILURE);
            }

    }
    break;
#endif


//______________________________________________________________________________

    case DPF_ABOUT:
        /*
        **  about
        **  Rotation center for current ligand,
        */
        nfields = sscanf( line, "%*s " FDFMT3, &lig_center[X], &lig_center[Y], &lig_center[Z]);
	if(nfields!=3) syntaxstop("ABOUT");
        pr( logFile, "Small molecule center of rotation =\t" );
        pr( logFile, "(%+.3f, %+.3f, %+.3f)\n\n", lig_center[X], lig_center[Y], lig_center[Z]);
        B_found_about_keyword = TRUE; //set false by 'move' true by 'about'
        B_found_tran0_keyword = FALSE;
        if ( nmol == 0 || true_ligand_atoms==0) {
            pr( logFile, "Must specify a ligand PDBQT file, using the \"move\" command.\n");

	    }
            /*
            **  Zero-out on central point...
            */
            maxrad = 0;
            for ( i=0; i<true_ligand_atoms; i++ ) { /*new, gmm, 6-23-1998*/
                r2sum=0.0;
                for (xyz = 0;  xyz < SPACE;  xyz++) {
		    Real c;
                    c = crdorig[i][xyz] - lig_center[xyz];
                    crdpdb[i][xyz] = c;
                    crd[i][xyz] = c;
                    r2sum += c*c;
                } /* xyz */
                maxrad = max(maxrad,sqrt(r2sum));
            } /* i */
            if (outlev >= LOGLIGREAD && true_ligand_atoms>0) {
                pr( logFile, "Furthest true ligand atom from \"about\" center is %.3f Angstroms (maxrad).\n\n",maxrad);
            }
        break;

/*____________________________________________________________________________*/

    case DPF_REORIENT:
        /*
         *  reorient random
         *      # applies a random rotation to the input ligand
         * -OR-
         *  reorient standard
         *      # moves the ligand such that
         *      # the first three atoms lie parallel to the xy-plane, and
         *      # the first two atoms lie parallel to the x-axis
         * -OR-
         *  reorient <axis-x> <axis-y> <axis-z> <angle>
         *      # applies the specified rotation to the input ligand
	 *
	 * this modifies crdpdb but does not touch crdorig
         */
        get1arg( line, "%*s %s", param[0], "REORIENT" );
        { // Parse the reorient command
            if (streq(param[0],"random")) {
                // reorient random
                B_reorient_random = TRUE; // create a new random orientation before docking
                q_reorient = randomQuat();

            } else if (streq(param[0],"standard")) {
                { // reorient standard
                B_reorient_random = FALSE; // do not create a new random orientation before docking

                if (true_ligand_atoms >= 3 ) {
                    // Move the ligand such that
                    // the first three atoms lie parallel to the xy-plane, and
                    // the first two atoms lie parallel to the x-axis
                    Vector vec_01,     // vector between ligand atoms 0 and 1
                           vec_12,     // vector between ligand atoms 1 and 2
                           vec_normal, // vector perpendicular to plane of vec_01 and vec_12
                           vec_x_axis, // vector along the X-axis
                           vec_z_axis, // vector along the Z-axis
                           vec_reorient_axis; // vector describing the axis about which to reorient
                    // Set the X and Z axes:
                    vec_x_axis[X] = 1.;
                    vec_x_axis[Y] = 0.;
                    vec_x_axis[Z] = 0.;
                    vec_z_axis[X] = 0.;
                    vec_z_axis[Y] = 0.;
                    for (xyz = 0;  xyz < SPACE;  xyz++) {
                        vec_01[xyz] = (double)( crdpdb[1][xyz] - crdpdb[0][xyz] );
                        vec_12[xyz] = (double)( crdpdb[2][xyz] - crdpdb[1][xyz] );
                    }
                    // Compute the normal to vec_01 and vec_12
                    Cross_product( vec_normal, vec_01, vec_12 );
                    Print_vector( logFile, "vec_01", vec_01 );
                    Print_vector( logFile, "vec_12", vec_12 );
                    Print_vector( logFile, "vec_normal", vec_normal );
                    Print_vector( logFile, "vec_z_axis", vec_z_axis );
                    // Compute the angle between vec_01 and vec_12
                    double angle_012 = 0.;
                    angle_012 = Angle_between( vec_01, vec_12 );
                    pr( logFile, "Angle between vectors 01 and 12 = %.2f degrees\n", RadiansToDegrees( angle_012 ) );
                    if ( ( fabs(angle_012) < APPROX_ZERO ) || ( ( fabs(angle_012) > (PI - APPROX_ZERO) ) && ( fabs(angle_012) < (PI + APPROX_ZERO) ) ) ) {
                        // angle is too small or "too linear" to align the molecule into the xy-plane
                        pr( logFile, "%s:  WARNING!  The angle between the first three atoms is not suitable (%6.3f degrees) to align them with the xy-plane.\n", programname, RadiansToDegrees( angle_012 ) );
                    } else {
                        // Calculate angle between vec_normal and the z-axis
                        double angle_n1z = 0.;  // Angle between vec_normal and the z-axis
                        angle_n1z = Angle_between( vec_normal, vec_z_axis );
                        pr( logFile, "Angle between vec_normal and vec_z_axis = %.2f degrees\n", RadiansToDegrees( angle_n1z ) );
                        //
                        // We need to rotate the molecule about the normal to vec_normal and vec_z_axis
                        Cross_product( vec_reorient_axis, vec_normal, vec_z_axis );
                        // Set the rotation axis for reorientation
                        // Set the angle for reorientation of the first 3 atoms
                        // into the xy-plane
                        q_reorient = raaDoubleToQuat(vec_reorient_axis, -angle_n1z);

                        // Rotate ligand into the xy-plane...
                        // qtransform( origin, q_reorient, crdpdb, true_ligand_atoms );
                        qtransform( origin, q_reorient, crdpdb, true_ligand_atoms );

                        // Compute the updated vec_01, the vector between atom 0 and atom 1,
                        // since the preceding "qtransform" changed the coordinates.
                        for (xyz = 0;  xyz < SPACE;  xyz++) {
                            vec_01[xyz] = (double)( crdpdb[1][xyz] - crdpdb[0][xyz] );
                        }
                        //
                        // Compute the angle between vec_01 and the x-axis:
                        double angle_01x = 0.;
                        angle_01x = Angle_between( vec_01, vec_x_axis );
                        //
                        pr( logFile, "Angle between vector 01 and the x-axis = %.2f degrees\n", RadiansToDegrees( angle_01x ) );
                        //
                        // The rotation axis to rotate the first two atoms, 0 and 1,
                        // to be parallel to the x-axis, will be
                        // perpendicular to the xy-plane, i.e. the z-axis,
                        // since the molecule's first 3 atoms are now in the xy-plane.
                        // Set the rotation angle:
                        // Build the quaternion from the axis-angle rotation values:
                        q_reorient = raaDoubleToQuat(vec_z_axis, angle_01x);
                    } // angle_012 is appropriate to align into xy-plane

                } else {
                    prStr( error_message, "%s: ERROR! Insufficient atoms in the ligand.  There must be at least three atoms in the ligand to use this command.\n", programname );
                    stop( error_message ); // exits
                }
                } // reorient standard
            } else {
                { // reorient <nx> <ny> <nz> <angle>
		    AxisAngle aa;
                    B_reorient_random = FALSE; // do not create a new random orientation before docking

                    // Read the specified initial orientation for the ligand
                    nfields = sscanf( line,"%*s %lf %lf %lf %lf",
		     &aa.nx, &aa.ny, &aa.nz, &aa.ang );
                    if ( nfields == 4 ) {
                        // Normalise the vector defining the axis of rotation:
                        // Make sure angle is in radians, and ranges from -PI to PI
                        aa.ang = DegreesToRadians( aa.ang ); // convert from degrees to radians
                        aa.ang = ModRad( aa.ang ); // wrap to range (0, 2*PI) using modulo 2*PI
                        aa.ang = WrpRad( aa.ang ); // wrap to range (-PI, PI)
                        pr( logFile, "After normalising the vector, and converting the angle to radians, the axis-angle rotation becomes ((%.3f, %.3f, %.3f), %.2f radians)\n",
                                aa.nx, aa.ny, aa.ny, aa.ang);
                        // Convert the rotation-about-axis components (nx,ny,nz,ang)
                        // to a rotation-quaternion (x,y,z,w):
                        q_reorient = AxisAngleToQuat(aa);
                    } else {
                        prStr( error_message, "%s: ERROR! Please specify the vector x,y,z and rotation angle (degrees) using four real numbers.\n", programname );
                        stop( error_message );
                    }
                } // reorient <nx> <ny> <nz> <angle>
            } // endif
        } // end parsing reorient command line

        reorient( logFile, true_ligand_atoms, atomstuff, crdpdb, charge, type,
                  parameterArray, q_reorient, origin, ntor, tlist, vt, &ligand,
		  debug, outlev);

        break;


//______________________________________________________________________________

    case DPF_TRAN0:
        /*
        **  tran0
        **  Initial_translation,
        */
        get1arg( line, "%*s %s", param[0],  "TRAN0");
        if (streq(param[0],"random")) {
            B_RandomTran0 = TRUE;
            ligand.S.T.x = sInit.T.x = random_range( info->lo[X], info->hi[X] );
            ligand.S.T.y = sInit.T.y = random_range( info->lo[Y], info->hi[Y] );
            ligand.S.T.z = sInit.T.z = random_range( info->lo[Z], info->hi[Z] );
        } else {
            B_RandomTran0 = FALSE;
            nfields = sscanf( line,"%*s %lf %lf %lf", &(sInit.T.x), &(sInit.T.y), &(sInit.T.z));
	    if(nfields!=3) stop("syntax error in TRAN0 X Y Z line");
            ligand.S.T.x = sInit.T.x;
            ligand.S.T.y = sInit.T.y;
            ligand.S.T.z = sInit.T.z;
        }
        B_found_tran0_keyword = TRUE;
        if (outlev >= LOGBASIC) {
            pr( logFile, "Initial translation =\t\t\t(%.3f, %.3f, %.3f) Angstroms\n", sInit.T.x, sInit.T.y, sInit.T.z );
        }
        break;

//______________________________________________________________________________

    case DPF_QUAT0:
    case DPF_AXISANGLE0:
    case DPF_QUATERNION0:
        /*
         * Handles both axisangle0 and quaternion0
         *
         *  axisangle0 1. 0. 0. 0.
         *  axisangle0 random
         *  ( quat0 <--- deprecated )
         *  Initial_quaternion, specified as an axis and angle
         *
         *  quaternion0 0. 0. 0. 1.
         *  quaternion0 random
         *  Initial_quaternion, specified as the four components (qx, qy, qz, qw)
         */
        {
        // Local Block...
        double a, b, c, d;
        get1arg( line, "%*s %s", param[0], "QUATERNION0 or AXISANGLE0");
        if (streq(param[0],"random")) {
            // Make a random initial quaternion,
            // and set the boolean B_RandomQuat0 to true,
            // so we can generate random quaternions in population-based methods.
            B_RandomQuat0 = TRUE;
            sInit.Q = randomQuat();
            if (outlev >= LOGBASIC) {
                pr( logFile, "Each run will begin with a new, random initial orientation.\n");
            }
        } else {
            // Read in the user-defined axis-angle values for the initial quaternion
            // and set the boolean B_RandomQuat0 to false,
            B_RandomQuat0 = FALSE;
            nfields = sscanf( line, "%*s %lf %lf %lf %lf", &a, &b, &c, &d);
	    if(nfields!=4) stop("syntax error in AXISANGLE0 or QUATERNION0 values ");
            sInit.Q = (dpf_keyword == DPF_QUATERNION0) ?
                      quatComponentsToQuat(a,b,c,d) :
                      axisDegreeToQuat(a,b,c,d);
        }

        ligand.S.Q = sInit.Q;
	// LOGTODO fix logic in the next mess MP
        if (outlev >= LOGBASIC ) {
            if (dpf_keyword == DPF_QUATERNION0) {
                pr( logFile, "Initial quaternion,  (x,y,z,w) =\t( %.3f, %.3f, %.3f, %.3f ),\n", sInit.Q.x, sInit.Q.y, sInit.Q.z, sInit.Q.w);
            } else {
                if (dpf_keyword == DPF_QUAT0 && !B_RandomQuat0)  {
                    pr( logFile, "WARNING quat0 command is obsolete. Now use quaternion0 or axisangle0 instead\n");
                }
                if (!B_RandomQuat0) {
                    pr( logFile, "Initial axis-angle,  (nx,ny,nz,ang) =\t( %.3f, %.3f, %.3f, %.1f deg ),\n", a, b, c, d );
                }
                pr( logFile, "Initial quaternion,  (x,y,z,w) =\t( %.3f, %.3f, %.3f, %.3f ),\n", sInit.Q.x, sInit.Q.y, sInit.Q.z, sInit.Q.w);
            }
#ifdef DEBUG
            pr( logFile, "Initial Quaternion sInit.Q:\n\n");
            printQuat( logFile, sInit.Q );
            pr( logFile, "Initial Quaternion ligand.S.Q:\n\n");
            printQuat( logFile, ligand.S.Q );
#endif
        }
        } // end Local Block
        break;

//______________________________________________________________________________

    case DPF_NDIHE:
        /*
        **  ndihe
        **  Formerly, number of dihedral angles to be specified by "dihe0"
        */
            if (outlev >= LOGMIN) {
                pr( logFile, "%s: WARNING!  The \"ndihe\" command is no longer supported.  The number of torsions in the PDBQT file(s) is the number that will be used (i.e. %d)\n", programname, ntor);
            }
        break;

//______________________________________________________________________________

    case DPF_DIHE0:
        /*
        **  dihe0
        **  Initial dihedral angles, input in degrees,
        */
        get1arg( line, "%*s %s", param[0], "DIHE0");
        if (streq(param[0],"random")) {
            B_RandomDihe0 = TRUE;
            sInit.ntor = nval = ntor;
            for ( i=0; i<nval; i++ ) {
                sInit.tor[i] = random_range( -180.0, 180.0 );
            }
        } else {
            B_RandomDihe0 = FALSE;
            nfields = sscanf( line, torfmt, TOR_ARG_LIST );
            if (nfields == 0) {
                stop( "could not read any torsions in DIHE0 line" );
            } else if (nfields == EOF) {
                stop( "End of file encountered while reading DIHE0 line");
            } else if (nfields < ntor) {
                pr( logFile, "Only %d initial torsion angles were detected on input DIHE0 line.\n",nfields);
                pr( logFile, "The number of torsions detected in the PDBQT files was %d torsions.\n", ntor);
	        stop("torsion count mismatch");
            } else {
                if (outlev >= LOGBASIC) {
                    pr( logFile, "%d initial torsion angles were detected on input line.\n", nfields );
                }
            }
            nval = nfields;
        }
        if (nval != ntor) {
            pr( logFile, "%s: WARNING!  The number of torsions specified (%d) does not match the number found in the PDBQT file (i.e. %d)\n", programname, nval, ntor);
        }
        for ( i=0; i<nval; i++ ) {
            if (outlev > LOGFORADT) 
                pr( logFile, "\tInitial torsion %2d = %7.2f deg\n", (i+1), sInit.tor[i] ); /* sInit.tor is in degrees */
                /* Convert sInit.tor[i] into radians */
            ligand.S.tor[i] = sInit.tor[i] = DegreesToRadians( sInit.tor[i] ); /* sInit.tor is now in radians  Added:05-01-95 */
        }
        break;

//______________________________________________________________________________

    case DPF_TSTEP:
        /*
        **  tstep
        **  Simulated annealing Translation_step,
        */
        nfields = sscanf( line, "%*s " FDFMT2, &trnStep0, &trnStepFinal );
        if (nfields == 0) {
            stop( " Could not read any arguments in TSTEP line" );
        } else if (nfields == EOF) {
            stop(  "End of file encountered in TSTEP line");
        } else if (nfields > 0) {
            pr( logFile, "Initial simanneal cycle, maximum translation step = +/- %-.1f Angstroms\n", trnStep0);
        }
        if (nfields == 2) {
            B_CalcTrnRF = TRUE;
            if (outlev >= LOGBASIC) {
                pr( logFile, "Final cycle,   maximum translation step = +/- %-.1f Angstroms\n", trnStepFinal);
                pr( logFile, "Reduction factor will be calculated when number of cycles has been read in.\n");
            }
        }
        break;

//______________________________________________________________________________

    case DPF_QSTEP:
        /*
        **  qstep
        **  Simulated annealing Quaternion_step,
        */
        nfields = sscanf( line, "%*s " FDFMT2, &qtwStep0, &qtwStepFinal );
        if (nfields == 0) {
            stop("could not read any arguments in QSTEP line" );
        } else if (nfields == EOF) {
            stop("End of file encountered in QSTEP line");
        } else if (nfields > 0) {
            if (outlev >= LOGBASIC) {
                pr( logFile, "Initial simanneal cycle, maximum quaternion angle step = +/- %-.1f deg\n", qtwStep0);
            }
            /* convert to radians */
            qtwStep0 = DegreesToRadians( qtwStep0 );
        }
        if (nfields == 2) {
            B_CalcQtwRF = TRUE;
            if (outlev >= LOGBASIC) {
                pr( logFile, "Final cycle,   maximum quaternion angle step = +/- %-.1f deg\n", qtwStepFinal);
                pr( logFile, "Reduction factor will be calculated when number of cycles has been read in.\n");
            }
            /* convert to radians */
            qtwStepFinal = DegreesToRadians( qtwStepFinal );
        }
        break;

//______________________________________________________________________________

    case DPF_DSTEP:
        /*
        **  dstep
        **  Simulated annealing Torsion_step,
        */
        nfields = sscanf( line, "%*s " FDFMT2, &torStep0, &torStepFinal );
        if (nfields == 0) {
            stop( "Could not read any arguments in DSTEP line" );
        } else if (nfields == EOF) {
            stop( "End of file encountered in DSTEP line");
        } else if (nfields > 0) {
            if (outlev >= LOGBASIC) {
                pr( logFile, "Initial simanneal cycle, maximum torsion angle step = +/- %-.1f deg\n", torStep0);
            }
            /* convert to radians */
            torStep0 = DegreesToRadians( torStep0 );
        }
        if (nfields == 2) {
            B_CalcTorRF = TRUE;
            if (outlev >= LOGBASIC) {
                pr( logFile, "Final simanneal cycle,   maximum torsion angle step = +/- %-.1f deg\n", torStepFinal);
                pr( logFile, "Reduction factor will be calculated when number of simanneal cycles has been read in.\n");
            }
            /* convert to radians */
            torStepFinal = DegreesToRadians( torStepFinal );
        }
        break;

//______________________________________________________________________________

    case DPF_TRNRF:
        /*
        **  trnrf
        **  Translation reduction factor,
        */
        get1arg( line, "%*s " FDFMT, &trnFac, "TRNRF" );
        if (outlev >= LOGBASIC) {
            pr( logFile, "Reduction factor for simanneal translations =\t%-.3f /cycle\n", trnFac );
        }
        B_trnReduc = (trnFac != 1.);
        break;

//______________________________________________________________________________

    case DPF_QUARF:
        /*
        **  quarf
        **  Quaternion reduction factor,
        */
        get1arg( line, "%*s " FDFMT, &qtwFac, "QRARF" );
        if (outlev >= LOGBASIC) {
            pr( logFile, "Reduction factor for simanneal quaternion angle =\t%-.3f /cycle\n", qtwFac );
        }
        B_qtwReduc = (qtwFac != 1.);
        break;

//______________________________________________________________________________

    case DPF_DIHRF:
        /*
        **  dihrf
        **  Torsion reduction factor,
        */
        get1arg( line, "%*s " FDFMT, &torFac, "DIHRF" );
        if (outlev >= LOGBASIC) {
            pr( logFile, "Reduction factor for simanneal torsion angles =\t%-.3f /cycle\n", torFac );
        }
        B_torReduc = (torFac != 1.);
        break;

//______________________________________________________________________________

    case DPF_FLEX:
        /*
        **  flex
        **  Flexible side-chains, cannot translate:
        */
        nmol++;
        nres++;
        break;

//______________________________________________________________________________

    case DPF_INTNBP_REQM_EPS:
    case DPF_INTNBP_COEFFS:
        /*
        **  intnbp_r_eps
        **  Read internal energy parameters:
        **  Lennard-Jones and Hydrogen Bond Potentials,
        **  DPF_INTNBP_REQM_EPS: Using epsilon and r-equilibrium values...
        **  DPF_INTNBP_COEFFS: Using coefficients...
        */
	{ // block for locals
        Real epsij;
        Real Rij;
	int xA, xB;
        nfields = sscanf( line, "%*s " FDFMT2 " %d %d %s %s", &Rij, &epsij, &xA, &xB, param[0], param[1] );
	if(nfields!=6) stop("syntax error, not 6 values in INTNBP_R_EPS line");
        if ( dpf_keyword == DPF_INTNBP_REQM_EPS ) {
        /* check that the Rij is reasonable */
	/* SF ...but only if there are no G-atoms. */        /* SF RING CLOSURE */

	if ((Rij <= 2.0 ) && (epsij >= EPSIJ_MAX )) {    /* RING CLOSURE */
 	     (void) fprintf( logFile, "Ring closure distance potential found for atom type %s :\n    Equilibrium distance   = %.2f Angstroms \n    Equilibrium potential  = %.6f Kcal/mol\n    Pseudo-LJ coefficients = %d-%d \n\n", param[1] , Rij, epsij, xA, xB); /* SF RING CLOSURE */
			}   /* SF RING CLOSURE */
	else { /* SF RING CLOSURE */

	        if ((Rij < RIJ_MIN) || (Rij > RIJ_MAX)) {
        	    (void) fprintf( logFile,
	            "WARNING: pairwise distance, Rij, %.2f, is not a very reasonable value for the equilibrium separation of two atoms! (%.2f Angstroms <= Rij <= %.2f Angstroms)\n\n", Rij, RIJ_MIN, RIJ_MAX);
	            (void) fprintf( logFile, "Perhaps you meant to use \"intnbp_coeffs\" instead of \"intnbp_r_eps\"?\n\n");
	            /* GMM COMMENTED OUT FOR DAVE GOODSELL, MUTABLE ATOMS
	             * exit(EXIT_FAILURE); */
	     	     }
	        /* check that the epsij is reasonable */
	        if ((epsij < EPSIJ_MIN) || (epsij > EPSIJ_MAX)) {
	            (void) fprintf( logFile,
	            "WARNING: well-depth, epsilon_ij, %.2f, is not a very reasonable value for the equilibrium potential energy of two atoms! (%.2f kcal/mol <= epsilon_ij <= %.2f kcal/mol)\n\n", epsij, EPSIJ_MIN, EPSIJ_MAX);
	            (void) fprintf( logFile, "Perhaps you meant to use \"intnbp_coeffs\" instead of \"intnbp_r_eps\"? \n\n");
	            /* GMM COMMENTED OUT FOR DAVE GOODSELL, MUTABLE ATOMS
	             * exit(EXIT_FAILURE); */
	        }
	        } /* RING CLOSURE */

	     }

        /* Defend against division by zero... */
        if (xA != xB) {
            if ( dpf_keyword == DPF_INTNBP_REQM_EPS ) {
               // Calculate the coefficients from Rij and epsij
    	       double tmpconst = epsij / (Real)(xA - xB);
               cA = tmpconst * pow( (double)Rij, (double)xA ) * (Real)xB;
               cB = tmpconst * pow( (double)Rij, (double)xB ) * (Real)xA;
            } else {
               cA = Rij;
               cB = epsij;
            }

            int a[2]; /* atom types of this interaction pair */
	    Boole is_hbond = FALSE;  // MPique 2012 not implemented here
            for (int i=0;i<2;i++) {
                foundParameter = apm_find(param[i]);
                if ( NULL == foundParameter ) {
                    prStr( error_message,"%s: ERROR:  Unknown ligand atom type \"%s\"; add parameters for it to the parameter library first!\n", programname, param[i]);
                    stop(" unknown ligand atom type");
		    /* NOTREACHED */
                }
                else a[i] = foundParameter->map_index;
            }
            pr(logFile, "\nCalculating internal non-bonded interaction energies for docking calculation;\n");
            intnbtable( &B_havenbp, a[0], a[1], info, cA, cB, xA, xB, is_hbond, r_smooth, AD4, sigma, ad_energy_tables, BOUND_CALCULATION, logFile, outlev);
           pr(logFile, "\nCalculating internal non-bonded interaction energies for unbound conformation calculation;\n");
           intnbtable( &B_havenbp, a[0], a[1], info, cA_unbound, cB_unbound, xA_unbound, xB_unbound, is_hbond, r_smooth, AD4, sigma, unbound_energy_tables, UNBOUND_CALCULATION, logFile, outlev );
        } else {
            stop("ERROR: Exponents must be different, to avoid division by zero!\n\tAborting...\n");
            exit(EXIT_FAILURE);
        }
	} // block for locals
        break;

//______________________________________________________________________________


    case DPF_UNBOUND_INTNBP_COEFFS:
        /*
        **  unbound_intnbp_coeffs
        **  Read internal energy parameters for unbound extended state calculation:
        */
        nfields = sscanf( line, "%*s " FDFMT2 " %d %d", &cA_unbound, &cB_unbound, &xA_unbound, &xB_unbound );
	if(nfields!=2) stop("syntax error, not 2 values in UNBOUND_INTNBP_COEFFS line");

        pr(logFile, "\nSetting the internal non-bonded interaction energy parameters for the\nunbound docking calculation, E = %.1f / r^%d - %.1f / r^%d\n\n", cA_unbound, xA_unbound, cB_unbound, xB_unbound);
        break;

//______________________________________________________________________________

    case DPF_RT0:
        /*
        **  rt0
        **  Initial Temperature,
        */
        get1arg( line, "%*s " FDFMT, &RT0, "RT0" );
        if (RT0 <= 0.) {
	    stop("Negative or zero temperature in RT0 line");
        }
        if (outlev >= 0) {
            pr( logFile, "\n\t\tTEMPERATURE SCHEDULE INFORMATION\n" );
            pr( logFile, "\t\t________________________________\n\n" );
            pr( logFile, "               -1 -1                 -1 -1\n" );
            pr( logFile, "R = %5.3f J mol  K    = %5.3f cal mol  K  \n\n", RJ, Rcal );
            pr( logFile, "                                        -1\n" );
            pr( logFile, "Initial R*Temperature = %8.2f cal mol\n", RT0 );
            pr( logFile, "      (=> Temperature = %8.2f K\n", RT0/Rcal );
            pr( logFile, "                   or = %8.2f C)\n\n", RT0/Rcal - T0K );
        }
        break;

//______________________________________________________________________________

    case DPF_RTRF:
        /*
        **  rtrf
        **  Temperature reduction factor,
        */
        get1arg( line, "%*s " FDFMT, &RTFac, "RTRF");
        if (outlev >= LOGBASIC) {
            pr( logFile, "R*Temperature reduction factor = %8.2f\t/cycle\n", RTFac );
        }
        if (RTFac >= 1.) {
            stop("Cooling is impossible with a reduction\n\tfactor greater than or equal to 1.0!" );
        } else if (RTFac == 0.0 ) {
            stop("Cooling is impossible with a ZERO reduction factor!" );
        } else if (RTFac < 0.0 ) {
            stop("Cooling is impossible with a NEGATIVE reduction factor!" );
        }
        break;

//______________________________________________________________________________

    case DPF_RUNS:
        /*
        **  runs
        **  Number of docking runs: GA or LS or simanneal
	**  Note this need not be checked here against MAX_RUNS-nconf as DPF could
	**  modify it before triggering runs M Pique
        */
        get1arg( line, "%*s %d", &nruns, "RUNS" );
        if ( nruns > MAX_RUNS ) {
            prStr( error_message, "%s:  ERROR: %d runs were requested, but AutoDock is only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", programname, nruns, MAX_RUNS);
            stop( error_message );
        }
        pr( logFile, "Number of runs = %d run%s\n", nruns, pl(nruns));
        break;

//______________________________________________________________________________

    case DPF_CYCLES:
        /*
        **  cycles
        **  Number of constant temperature SA cycles,
        */
        get1arg( line, "%*s %d", &ncycles, "CYCLES" );
        if (ncycles < 0) stop("Negative number of cycles in CYCLES line");
	if (outlev >= LOGBASIC)  {
           pr( logFile, "Maximum number of cycles = %8d cycles\n", ncycles);
        }
        break;

//______________________________________________________________________________

    case DPF_ACCS:
        /*
        **  accs
        **  Maximum number of simanneal steps accepted,
        */
        get1arg( line, "%*s %d", &naccmax, "ACCS" );
        if (naccmax < 0) {
            stop("Negative number of accepted moves in ACCS line");
        }
        if (outlev >= LOGBASIC) {
            pr( logFile, "Maximum number accepted per cycle =\t\t%8d steps\n", naccmax);
        }
        break;

//______________________________________________________________________________

    case DPF_REJS:
        /*
        **  rejs
        **  Maximum number of simanneal steps rejected,
        */
        get1arg( line, "%*s %d", &nrejmax, "REJS" );
        if (nrejmax < 0) stop("Negative number of rejected moves in REJS line");
        if (outlev >= LOGBASIC) {
            pr( logFile, "Maximum number rejected per cycle =\t\t%8d steps\n", nrejmax);
        }
        break;

//______________________________________________________________________________

    case DPF_SELECT:
        /*
        **  select
        **  Select either minimum or last state from previous simanneal cycle,
        */
	{ // block for locals
        char selminpar = 'm';
        get1arg( line, "%*s %c", &selminpar, "SELECT" );
	switch(selminpar) {
	  case 'm': case 'M': B_selectmin = TRUE; break;
	  case 'l': case 'L': B_selectmin = FALSE; break;
	  default: stop("unrecognized option in 'select' : must be 'l' or 'm'");
	  }
	} // block for locals
        break;

//______________________________________________________________________________

    case DPF_RMSTOL:
        /*
        **  rmstol
        **  Cluster tolerance,
        */
        get1arg( line, "%*s " FDFMT, &clus_rms_tol, "RMSTOL");
        if (outlev >= LOGBASIC) {
            pr( logFile, "Maximum RMS tolerance for conformational cluster analysis = %.2f Angstroms\n", clus_rms_tol);
        }
        break;

//______________________________________________________________________________

    case DPF_RMSREF:
        /*
        **  rmsref
        **  RMS Reference Coordinates:
        */
        get1arg( line, "%*s %s", FN_rms_ref_crds, "RMSREF");
        if (outlev >= LOGBASIC) {
            pr( logFile, "RMS reference coordinates will taken from \"%s\"\n", FN_rms_ref_crds );
        }
        break;

//______________________________________________________________________________

    case DPF_RMSATOMS:
        /*
        **  rmsatoms ligand_only
        **  rmsatoms all
        **
        **  Set the atoms to compute the RMSD values for cluster analysis
        **  either "ligand_only" (the default) or "all" moving atoms (ligand + receptor)
        */
        nfields = sscanf( line, "%*s %s", rms_atoms_cmd);
        if (nfields != 1) {
            pr( logFile, "%s:  ERROR: please specify an argument (either \"ligand_only\" or \"all\").  By default, only the ligand atoms will be used for the cluster analysis.\n", programname );
	    stop("error in RMSATOMS line");
        } else {
            if ( streq( rms_atoms_cmd, "ligand_only")) {
                if (outlev >= LOGBASIC) {
                    pr( logFile, "RMS clustering will be performed on the ligand atoms only.\n" );
                }
                B_rms_atoms_ligand_only = TRUE;  // cluster on the ligand atoms only
            } else if ( streq( rms_atoms_cmd, "all")) {
                if (outlev >= LOGBASIC) {
                    pr( logFile, "RMS clustering will be performed on the moving atoms of the receptor plus all the ligand atoms.\n" );
                }
                B_rms_atoms_ligand_only = FALSE;  // cluster on the ligand atoms plus moving receptor atoms
            } else {
                if (outlev >= LOGBASIC) {
                    pr( logFile, "RMS clustering will be performed on the ligand atoms only.\n" );
                }
                B_rms_atoms_ligand_only = TRUE;  // cluster on the ligand atoms only
            }
        }
        break;

//______________________________________________________________________________

    case DPF_TRJFRQ:
        /*
        **  trjfrq
        **  Trajectory frequency,
        */
        get1arg( line, "%*s %d", &trj_freq, "TRJFRQ");
        B_write_trj = (trj_freq > 0);
        if (outlev >= LOGBASIC) {
            pr( logFile, UnderLine );
            pr( logFile, "\t\tTRAJECTORY INFORMATION\n" );
            pr( logFile, "\t\t______________________\n\n\n" );
        }
        if (B_write_trj) {
            if (outlev >= LOGBASIC) {
                pr( logFile, "Output frequency for simanneal trajectory frames =\tevery %d step%s\n", trj_freq, (trj_freq > 1)?"s.":"." );
            }
        } else {
            if (outlev >= LOGBASIC) {
                pr( logFile, "No trajectory of simanneal states will be written.\n\n" );
                pr( logFile, "Subsequent \"trjbeg\", \"trjend\", \"trjout\" and \"trjsel\" parameters will be ignored.\n\n" );
            }
        }
        break;

//______________________________________________________________________________

    case DPF_TRJBEG:
        /*
        **  trjbeg
        **  Trajectory begin cycle,
        */
        get1arg( line, "%*s %d", &trj_begin_cyc, "TRJBEG" );
        if (outlev >= LOGBASIC) {
            pr( logFile, "Begin outputting trajectory of states at cycle:\t%d\n", trj_begin_cyc );
        }
        if (trj_begin_cyc < 0) {
            trj_begin_cyc = 0;
        } else if (trj_begin_cyc > ncycles) {
            trj_begin_cyc = trj_end_cyc = ncycles;
        }
        --trj_begin_cyc;
        break;

//______________________________________________________________________________

    case DPF_TRJEND:
        /*
        **  trjend
        **  Trajectory end cycle,
        */
        get1arg( line, "%*s %d", &trj_end_cyc, "TRJEND" );
        if (outlev >= LOGBASIC) {
            pr( logFile, "Cease outputting trajectory of states at cycle:\t%d\n", trj_end_cyc );
        }
        if (trj_end_cyc > ncycles) {
            trj_end_cyc = ncycles;
        } else if (trj_end_cyc < 0) {
            trj_end_cyc = 1;
        }
        --trj_end_cyc;
        break;

//______________________________________________________________________________

    case DPF_TRJOUT:
        /*
        **  trjout
        **  Trajectory file,
        */
        get1arg( line, "%*s %s", FN_trj, "TRJOUT" );
        if (outlev >= LOGBASIC) {
            pr( logFile, "\nWrite trajectory of state variables to file: \"%s\"\n", FN_trj);
        }
        break;

//______________________________________________________________________________

    case DPF_TRJSEL:
        /*
        **  trjsel
        **  Trajectory select,
        */
	{ // block for locals
        char out_acc_rej = '?';
        get1arg( line, "%*s %c", &out_acc_rej, "TRJSEL" );
        B_acconly = (out_acc_rej == 'A');
        B_either  = (out_acc_rej == 'E');
	if(!(B_acconly || B_either)) {
            stop("Missing or unknown accepted/rejected TRJSEL output flag.\n" );
        }
	} // block for locals
        break;

//______________________________________________________________________________

    case DPF_EXTNRG:
        /*
        **  extnrg
        **  Wall Energy,
        */
        get1arg( line, "%*s " FDFMT, &WallEnergy, "EXTNRG" );
        if (outlev >= LOGBASIC) {
            pr( logFile, "External grid energy (beyond grid map walls) = %.2f\n\n", WallEnergy );
        }
        break;

//______________________________________________________________________________

    case DPF_CLUSTER:
        /*
        **  cluster
        **  Cluster mode,
        */
        get1arg( line, "%*s %s", FN_clus, "CLUSTER" );
        B_cluster_mode = TRUE;
        if (outlev >= LOGBASIC) {
            pr( logFile, "Cluster mode is now set.\n\n" );
        }
        clmode( num_atom_types, clus_rms_tol,
                hostnm, jobStart, tms_jobStart,
                B_write_all_clusmem, FN_clus, crdpdb, lig_center,
                B_symmetry_flag, B_unique_pair_flag, FN_rms_ref_crds,
                B_rms_heavy_atoms_only, h_index, outlev, logFile);
	// note : clmode exits program !
        break;

//______________________________________________________________________________

    case DPF_CLUSALL:
        /*
        ** write_all_clusmem
        ** Write all cluster members...
        */
        B_write_all_clusmem = TRUE;
        if (outlev >= LOGBASIC) {
            pr( logFile, "All members of each cluster will be written out after the clustering histogram.\n(This is instead of outputting just the lowest energy member in each.)\n\n" );
        }
        break;

//______________________________________________________________________________

    case DPF_RMSNOSYM:
        /*
        **  rmsnosym
        **  Calculate RMS values in the normal way,
        **  ignoring any atom-type equivalences...
        */
        B_symmetry_flag = FALSE;
        if (outlev >= LOGBASIC) {
            pr( logFile, "Symmetry will be ignored in RMS calculations.\n\n" );
        }
        break;

//______________________________________________________________________________

    case DPF_RMSMODE:
        /*
        **  rmsmode atype|unique_pair|heavy_atoms_only
        **  Calculate symmetrical RMS values 
        **  considering each atom to be paired at most one time (unique_pair)
        **  considering only non-hydrogen atoms (heavy_atoms)
        **  or used repeatedly (atype) (AD 4.2 default)
        */
        char rms_mode[LINE_LEN];
        get1arg( line, "%*s %s", rms_mode , "RMS_MODE");
        if (streq(rms_mode, "unique_pair")||streq(rms_mode, "uniquepair")) {
            B_unique_pair_flag = TRUE;
            if (outlev >= LOGBASIC) {
                pr( logFile, "Symmetry in RMS calculations will consider only unique atom pairs.\n\n" );
            }
        } else if (streq(rms_mode, "atype")) {
            B_unique_pair_flag = FALSE;
            if (outlev >= LOGBASIC) {
                pr( logFile, "Symmetry in RMS calculations will consider all atom pairs.\n\n" );
            }
        } else if (streq(rms_mode, "heavy_atoms_only")) {
            B_rms_heavy_atoms_only = TRUE;  // cluster on the ligand heavy atoms only, excluding hydrogens
            if (outlev >= LOGBASIC) {
               pr( logFile, "RMS calculations will consider only heavy atom pairs.\n\n" );
            }
        } else {
            pr( logFile, "%s:  ERROR:  Unrecognized rms mode type \"%s\" .\n",
                    programname, rms_mode );
            stop("");
        }
        if (B_rms_heavy_atoms_only && B_unique_pair_flag)
            if (outlev >= LOGBASIC) {
               pr( logFile, "RMS calculations will consider only unique pairs of heavy atoms.\n\n" );
            }
        break;

//______________________________________________________________________________

    case DPF_SCHEDGEOMETRIC:
        /*
        **  geometric_schedule
        **  Use a deprecated geometric temperature
        **  reduction schedule.  This was the default before 4.2.5
        */
        B_linear_schedule = FALSE;
        if (outlev >= LOGBASIC) {
            pr( logFile, "A geometric temperature reduction schedule will be used...\n\n" );
        }
        break;

//______________________________________________________________________________
//______________________________________________________________________________

    case DPF_SCHEDLIN:
        /*
        **  linear_schedule
        **  Use a linear (arithmetic) temperature
        **  reduction schedule.  This is necessary for
        **  more accurate entropy estimations...
	**  This is the default as of 4.2.5
        */
        B_linear_schedule = TRUE;
        if (outlev >= LOGBASIC) {
            pr( logFile, "A linear temperature reduction schedule will be used...\n\n" );
        }
        break;

//______________________________________________________________________________

    case DPF_WATCH:
        /*
        **  watch
        **  for watching a job's simanneal progress PDBQT file in AVS,
        */
        get1arg( line, "%*s %s", FN_watch, "WATCH");
        if (B_write_trj) {
            pr(logFile,"\nAutoDock will create the simanneal watch-file \"%s\", for real-time monitoring of runs.\n\n", FN_watch);
            pr(logFile,"\nThe watch-file will be updated every %d moves, in accordance with the trajectory parameters..\n\n", trj_freq);
            B_watch = TRUE;
        } else {
            pr(logFile,"\nYou must set \"trjfrq\" to be greater than zero. No watch-file will be created.\n\n");
            B_watch = FALSE;
        }
        break;

//______________________________________________________________________________

    case DPF_GAUSSTORCON:
    case DPF_HARDTORCON:
        /*
        ** "gausstorcon" Add Gaussian torsion contraints,
        ** "hardtorcon"  Add Hard torsion contraints,
        */
        nfields = sscanf( line, "%*s %d " FDFMT2, &I_tor, &F_torPref, &F_torHWdth);
	// I am not sure how many tokens are needed, so not guarding. M Pique 2010
        if (I_tor <= 0) {
            pr( logFile, "\nTorsion IDs less than 1 (%d) are not allowed!\n\n", I_tor);
	    stop("");
        } else if (I_tor > ntor) {
            pr( logFile, "\nRequested torsion ID (%d) is larger than the number of torsions found (%d)!\n\n", I_tor, ntor);
        } else { /* torsion-ID accepted */
            --I_tor;    /* Because humans start at 1, and C at 0... */
            if ( B_isTorConstrained[I_tor] == 0 ) {

                if (dpf_keyword ==  DPF_GAUSSTORCON) {
                    B_isGaussTorCon = TRUE;
                    B_isTorConstrained[I_tor] = 1;
                    /*
                    ** Initialize... Torsion Energy Profile...
                    ** Set energies at every torsion division
                    ** to the user-defined (maximum) barrier energy,
                    */
                    for (US_tD = 0;  US_tD < NTORDIVS;  US_tD++) {
                        US_torProfile[I_tor][US_tD] = US_torBarrier;
                    }
                } else {
                    /*
                    ** DPF_HARDTORCON
                    */
                    B_isTorConstrained[I_tor] = 2;
                }
            }
            if (dpf_keyword ==  DPF_GAUSSTORCON) {
                (void) strcpy( S_contype, " half-" );
            } else {
                (void) strcpy( S_contype, " " );
            }
                /*
            ** F_torPref ranges from -180.0 to +180.0 degrees...
            */
            F_torPref = WrpDeg(ModDeg(F_torPref));
            if (F_torHWdth < 0.) {
                pr(logFile,"\nI'm sorry, negative%swidths (%.1f) are not allowed. I will use the default (%.1f) instead.\n\n", S_contype, F_torHWdth, DEFHWDTH);
                F_torHWdth = DEFHWDTH;
            } else if (F_torHWdth > 90.) {
                pr(logFile,"\nI'm sorry, your requested%swidth (%.1f) is too large. I will use the default (%.1f) instead.\n\n", S_contype, F_torHWdth, DEFHWDTH);
                F_torHWdth = DEFHWDTH;
            }
            pr( logFile, "For torsion %d, Adding a constrained-torsion zone centered on %.1f degrees;\n%swidth = %.1f degrees.\n\n", 1+I_tor, F_torPref, S_contype, F_torHWdth);

            if (dpf_keyword == DPF_GAUSSTORCON) {
                /*
                ** Calculate the torsion energy profile;
                ** combine this with previous profile without
                ** losing any information.
                */
                for (F_A = F_A_from;  F_A <= F_A_to;  F_A += F_W) {
                    F_Aova = (F_A - F_torPref) / F_torHWdth;
                    US_energy = (unsigned short) (((Real)US_torBarrier) * (1.0 - exp(F_lnH * F_Aova*F_Aova)));
                    /*
                    ** if F_A(<-180.or>180), wrap to -180to180,
                    */
                    F_tor = WrpDeg(ModDeg(F_A));
                    /*
                    ** Convert from F_tor to US_tD
                    */
                    US_tD = (unsigned short) ((F_tor - F_hW + 180.)/F_W);
                    US_torProfile[I_tor][US_tD] = min(US_energy,US_torProfile[I_tor][US_tD]);
                }/* for F_A */
                /*
                ** Ensure that the lowest point(s) in the profile are
                ** zero...
                */
                US_min = TORBARMAX;
                for (US_tD = 0;  US_tD < NTORDIVS;  US_tD++) {
                    US_min = min(US_min,US_torProfile[I_tor][US_tD]);
                }
                for (US_tD = 0;  US_tD < NTORDIVS;  US_tD++) {
                    US_torProfile[I_tor][US_tD] -= US_min;
                }
            } else { /*DPF_HARDTORCON*/

                iCon = N_con[I_tor] + 1;
                if (iCon < MAX_TOR_CON) {
                    F_TorConRange[I_tor][N_con[I_tor]][LOWER] = F_torPref - 0.5* F_torHWdth;
                    F_TorConRange[I_tor][N_con[I_tor]][UPPER] = F_torPref + 0.5* F_torHWdth;
                    N_con[I_tor] = iCon;
                } else {
                    pr(logFile,"\n\n I'm sorry, you can only have %d (=MAX_TOR_CON) torsion constraints.\nIf you need more, change the \"#define MAX_TOR_CON\" line in \"constants.h\".\n\n",MAX_TOR_CON);
                }/* Still room to add another constraint. */
            } /*DPF_HARDTORCON*/
        }/* torsion-ID accepted */
        break;

//______________________________________________________________________________

    case DPF_BARRIER:
        /*
        **  barrier
        **  Define torsion-barrier energy...
        */
        get1arg( line, "%*s %d", &I_torBarrier, "BARRIER");
        US_torBarrier = (unsigned short)I_torBarrier;
        US_torBarrier = min(US_torBarrier, TORBARMAX);
        pr(logFile,"\nTorsion barrier energy is set to %uhd\n\n", US_torBarrier);
        break;

//______________________________________________________________________________

    case DPF_SHOWTORPEN:
        /*
        **  showtorpen
        **  Show torsion's penalty energy.
        */
        B_ShowTorE = TRUE;
        pr(logFile,"\nConstrained torsion penalty energies will be stored during docking, and output after each run\n\n");
        break;

//______________________________________________________________________________

    case DPF_E0MAX:
        /*
        **  e0max
        **  Set simanneal max initial energy, and, optionally, number of retries
        */
        nfields = sscanf( line, "%*s " FDFMT " %d", &e0max, &MaxRetries );
        if (nfields == 0) {
            stop("Could not read any arguments in E0MAX line" );
        } else if (nfields == EOF) {
            stop("End of file encountered in E0MAX line");
        } else if (nfields == 1) {
            pr(logFile,"Using the default maximum number of retries for simanneal initialization, %d retries.\n", MaxRetries);
        } else if (nfields == 2) {
            pr(logFile,"Using user-specified maximum number of retries for simanneal initialization, %d retries.\n", MaxRetries);
        }
        if (e0max < 0.) stop("e0max must be positive");
        pr(logFile,"If the simanneal initial energy is greater than e0max, %.3f,\nthen a new, random initial state will be created.\n\n",e0max);
        break;

//______________________________________________________________________________

    case DPF_SIMANNEAL:
        /*
        ** simanneal
        */
	    /* optional argument is alternate way to specify nruns : */
    	    nfields = sscanf( line, "%*s %d",&nruns);

            if ( nruns+nconf > MAX_RUNS ) {
                prStr( error_message, "%s:  ERROR: %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", 
		programname, nruns+nconf, MAX_RUNS);
                stop( error_message );
		}
            pr( logFile, "Number of simanneal runs  = %d nruns\n\n", nruns);
            pr( logFile, "Maximum number of cycles per run = %d cycles\n\n", ncycles);
	    if (B_linear_schedule) {
	       RTreduc = RT0 / ncycles;
	       if (outlev >= LOGBASIC) {
		  pr( logFile, "A linear temperature reduction schedule was requested...\n" );
		  pr( logFile, "Annealing temperature will be reduced by %.3f cal mol per cycle.\n", RTreduc );
		  }
		}

        /*
        ** Calculate reduction factor based on initial and final step values,
        ** and number of cycles...
        */
            if (B_CalcTrnRF) {
                trnFac = RedFac(trnStep0, trnStepFinal, ncycles-1);
                pr( logFile, "Calculated reduction factor for simanneal translations     = %-.3f /cycle\n", trnFac);
            }
            B_trnReduc = (trnFac != 1.);
            if (B_CalcQtwRF) {
                qtwFac = RedFac(qtwStep0, qtwStepFinal, ncycles-1);
                pr( logFile, "Calculated reduction factor for simanneal quaternion angle = %-.3f /cycle\n", qtwFac );
            }
            B_qtwReduc = (qtwFac != 1.);
            if (B_CalcTorRF) {
                torFac    = RedFac(torStep0, torStepFinal, ncycles-1);
                pr( logFile, "Calculated reduction factor for simanneal torsion angles   = %-.3f /cycle\n", torFac );
            }
            B_torReduc = (torFac != 1.);
            B_tempChange = ( RTFac != 1.0 );

            if(outlev >= LOGBASIC) {
		if ( B_selectmin ) {
			pr( logFile, "%s will begin each new cycle\nwith the state of minimum energy from the previous annealing cycle.\n", programname);
		  } else {
			pr( logFile, "%s will begin each new cycle\nwith the last state from the previous annealing cycle.\n", programname);
		}
	    }
            pr(logFile, "\n");
            /*
            ** Number of ligands read in...
            */
            if (nlig > 0) {
                pr( logFile, "Total number of ligands read in by the DPF \"move\" command = %d\n", nlig );
            }
            if (nres > 0) {
                pr( logFile, "Total number of residues read in by the DPF \"flex\" command = %d\n", nres );
            }
            if (outlev >= LOGBASIC) {
                   pr( logFile, "Maximum possible number of steps per cycle = %d steps\n", naccmax+nrejmax);
            }

	    if (outlev >= LOGBASIC) {
                if (B_acconly) pr( logFile, "Output *accepted* states only.\n" );
                else if (B_either) pr( logFile, "Output *either* accepted or rejected states.\n" );
            } 

	    // set lig_center if not already set, use to center "crdpdb" ligand
	    center_ligand(crdorig, !B_found_about_keyword, natom, true_ligand_atoms,
	      tlist, ntor, crdpdb, lig_center, &sInit.T, &ligand.S.T,
	      outlev>=LOGBASIC, outlev, logFile);

            if (B_havenbp && outlev>=LOGNBINTEV)  nbe( info, ad_energy_tables, num_atom_types, outlev, logFile );
            if (B_cluster_mode) {
                clmode( num_atom_types, clus_rms_tol,
                        hostnm, jobStart, tms_jobStart,
                        B_write_all_clusmem, FN_clus, crdpdb, lig_center,
                        B_symmetry_flag, B_unique_pair_flag, FN_rms_ref_crds,
                        B_rms_heavy_atoms_only, h_index, outlev, logFile);
            }
            for (j = nconf; j < MAX_RUNS; j++) {
                econf[j] = torsFreeEnergy;
            }
            if (ad4_unbound_model==Unbound_Default) ad4_unbound_model = Unbound_Same_As_Bound;
            pr(logFile, "Unbound model to be used is %s.\n", report_parameter_library());
            /* ___________________________________________________________________
            **
            ** Begin the automated docking simulation,
	    **  using simulated annealing
            ** ___________________________________________________________________
            */
            simanneal( &nconf, Nnb, Nnb_array, &group_energy, true_ligand_atoms,
	    WallEnergy, atomstuff, charge, abs_charge, qsp_abs_charge, B_calcIntElec,
                        crd, crdpdb, dock_param_fn,
                        ad_energy_tables,
                        econf, B_either,
                        peratomE,
                        ncycles, nruns, runseed, jobStart,
                        map,
                        naccmax, natom, nonbondlist, nrejmax, ntor, 
                        sInit, sHist,   qtwFac, B_qtwReduc, qtwStep0,
                        B_selectmin, FN_ligand, lig_center, RT0, B_tempChange, RTFac,
                        tms_jobStart, tlist, torFac, B_torReduc, torStep0,
                        FN_trj, trj_end_cyc, trj_begin_cyc, trj_freq, trnFac,
                        B_trnReduc, trnStep0, type, vt, B_write_trj,
                        B_constrain_dist, atomC1, atomC2, sqlower, squpper,
                        B_linear_schedule, RTreduc,
                        B_watch, FN_watch,
                        B_isGaussTorCon, US_torProfile, B_isTorConstrained,
                        B_ShowTorE, US_TorE, F_TorConRange, N_con,
                        B_RandomTran0, B_RandomQuat0, B_RandomDihe0,
                        e0max, torsFreeEnergy, MaxRetries, ligand_is_inhibitor,
                        ignore_inter,
                        B_include_1_4_interactions, scale_1_4, scale_eintermol,
                        parameterArray, unbound_internal_FE,
                        info, B_use_non_bond_cutoff,
                        B_have_flexible_residues,
                        PDBQT_record,
                        ad4_unbound_model,
			outlev,
			logFile
                        );

        break;

//______________________________________________________________________________

    case DPF_SET_GA:

      if (GlobalSearchMethod != NULL) {
          pr(logFile, "Deleting the previous settings for the Genetic Algorithm.\n");
          (void) fflush(logFile);
          delete GlobalSearchMethod;
          GlobalSearchMethod = NULL;
      }

      if (debug > 0) {
        if(output_pop_stats.everyNgens>0) pr( logFile, "\n\tOutput population statistics every %u generations.\n", output_pop_stats.everyNgens );
        else pr( logFile, "\n\tNever output generation-based population statistics.\n");
      }
      GlobalSearchMethod = new Genetic_Algorithm(e_mode, s_mode, c_mode, w_mode, elitism, c_rate, m_rate, localsearch_freq,
                                                 window_size, num_evals, num_generations, output_pop_stats);
      // note: the trn/qtw/torStep0 values appear unused beyond gs.h 
      // I do not know whether low, high are used but changing them has no 
      //    obvious effect
      // According to gs.cc, alpha and beta are unused
      //   - M Pique  April 2012
      ((Genetic_Algorithm *)GlobalSearchMethod)->mutation_values( low, high, alpha, beta, trnStep0, qtwStep0, torStep0  );
      ((Genetic_Algorithm *)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor, outlev, logFile);

      if (s_mode==LinearRanking){
        (void)((Genetic_Algorithm *)GlobalSearchMethod)->set_linear_ranking_selection_probability_ratio(linear_ranking_selection_probability_ratio);
          pr( logFile, "\n\tSet linear_ranking_selection_probability_ratio to %f.\n", linear_ranking_selection_probability_ratio );
      }
      

      break;
//______________________________________________________________________________

    case DPF_SET_SW1:

      if (LocalSearchMethod != NULL) {
          pr(logFile, "Deleting the previous settings for the local search Solis-Wets algorithm (SW1 object).\n");
          delete LocalSearchMethod;
          LocalSearchMethod = NULL;
      }

      pr(logFile, "Creating a new Local Search object using the Solis-Wets algorithm (SW1) with the current settings.\n\n");
      LocalSearchMethod = new Solis_Wets1(7+sInit.ntor, max_its, max_succ, max_fail, rho, lb_rho, 2.0, 0.5);

      break;
//______________________________________________________________________________

    case DPF_SET_PSW1:

      if (LocalSearchMethod != NULL) {
          pr(logFile, "Deleting the previous settings for the local search pseudo-Solis-Wets algorithm (pSW1 object).\n");
          delete LocalSearchMethod;
          LocalSearchMethod = NULL;
      }

      pr(logFile, "Creating a new Local Search object using the pseudo-Solis-Wets algorithm (pSW1) with the current settings.\n\n");

      //  Allocate space for the variable rho's
      rho_ptr = new Real[7+sInit.ntor];
      lb_rho_ptr = new Real[7+sInit.ntor];

      //  Initialize the rho's corresponding to the translation
      //  0,1,2   x,y,z
      //  3,4,5,6 qx,qy,qz,qw
      //  7,...   tor1
//These scale values can be changed in the dpf, officially unsupported 4/2009
//Real psw_trans_scale = 1.0;//1 angstrom
//Real psw_rot_scale = 0.05;  //about 3 degrees, we think
//Real psw_tors_scale = 0.1; //about 6 degrees

      for (j=0; j<3; j++) {
         // j=0,1,2
         rho_ptr[j] = rho * psw_trans_scale;// formerly trnStep0;
         lb_rho_ptr[j] = lb_rho * psw_trans_scale; //once trnStepFinal;
      }

      //  Initialize the rho's corresponding to the quaterion
      for (; j<7; j++) {
         // j=3,4,5,6
         rho_ptr[j] = rho * psw_rot_scale;// formerly qtwStep0;
         lb_rho_ptr[j] = lb_rho * psw_rot_scale; //once qtwStepFinal;
      }

      //  Initialize the rho's corresponding to the torsions
      for (; j<7+sInit.ntor; j++) {
         // j=7,...
         rho_ptr[j] = rho * psw_tors_scale;// formerly torStep0;
         lb_rho_ptr[j] = lb_rho * psw_tors_scale;//formerly torStepFinal;
      }

      LocalSearchMethod = new Pseudo_Solis_Wets1(7+sInit.ntor, max_its, max_succ, max_fail, 2.0, 0.5, rho_ptr, lb_rho_ptr);

      break;

//______________________________________________________________________________

    case DPF_SET_PATTERN:

      if (LocalSearchMethod != NULL) {
          pr(logFile, "Deleting the previous settings for the local search Pattern Search algorithm (PS object).\n");
          delete LocalSearchMethod;
          LocalSearchMethod = NULL;
      }


      pr(logFile, "Creating a new Local Search object using the Pattern Search algorithm (PS) with the current settings.\n\n");
      LocalSearchMethod = new Pattern_Search(7+sInit.ntor, max_succ, rho, lb_rho, 2.0, 0.5, localsearch_freq);

      break;

//______________________________________________________________________________

    case DPF_GS:
    case DPF_GALS:
        (void) fflush( logFile );
        /*
        ** Genetic Algorithm-Local search,  a.k.a. Lamarckian Genetic Algorithm
        */
            nfields= sscanf(line, "%*s %d",&nruns); // optional way to specify nruns
            if ( nruns+nconf > MAX_RUNS ) {
                prStr( error_message, "%s:  ERROR: %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", programname, nruns+nconf, MAX_RUNS);
                stop( error_message );
            }  
	    if (GlobalSearchMethod==NULL) {
               prStr(error_message, "%s:  ERROR:  You must use \"set_ga\" to allocate a Global Optimization method.\n", programname);
                 stop(error_message);
             }
	    if (dpf_keyword==DPF_GALS && LocalSearchMethod==NULL) {
               prStr(error_message, "%s:  ERROR:  You must use \"set_psw1\" to allocate a Local Optimization method.\n", programname);
                 stop(error_message);
             }
            exit_if_missing_elecmap_desolvmap_about("gals");

	    // set lig_center if not already set, use to center "crdpdb" ligand
	    center_ligand(crdorig, !B_found_about_keyword, natom, true_ligand_atoms,
	      tlist, ntor, crdpdb, lig_center, &sInit.Center, &ligand.S.Center,
	      outlev>=LOGBASIC, outlev, logFile);

	    // save centered crdpdb coords as crd (not sure is needed - MP 2012
	    for(int a=0;a<natom;a++) for(xyz=0;xyz<SPACE;xyz++)  
	      crd[a][xyz]=crdpdb[a][xyz];

	    // set 'tran0' vector to same as 'about' if not specified (2011-09)
	    if ( ! B_found_tran0_keyword ) {
		    ligand.S.T = sInit.T = sInit.Center;
		    pr( logFile, 
		    "Setting 'tran0' value to same as 'about' value: %.3f %.3f %.3f\n",
		    ligand.S.T.x, ligand.S.T.y, ligand.S.T.z);
		    B_found_tran0_keyword = TRUE;
	    }

            pr( logFile, "Number of requested %s dockings = %d run%s\n", GlobalSearchMethod->shortname(), nruns, pl(nruns));
            if (ad4_unbound_model==Unbound_Default) ad4_unbound_model = Unbound_Same_As_Bound;
            pr(logFile, "Unbound model to be used is %s.\n", report_parameter_library());

#ifdef DEBUG
            pr( logFile, "\nAbout to call evaluate.setup(), sInit:\n\n");
            printState( logFile, sInit, 2 );
#endif


            evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom,
                            info, map, peratomE, nonbondlist, ad_energy_tables, Nnb,
                            Nnb_array, &group_energy,
			    B_calcIntElec, B_isGaussTorCon, B_isTorConstrained,
                            B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, sInit, ligand,
                            ignore_inter,
                            B_include_1_4_interactions, scale_1_4, scale_eintermol,
                            unbound_internal_FE, B_use_non_bond_cutoff, B_have_flexible_residues,
			    true_ligand_atoms, outlev, logFile);

            evaluate.compute_intermol_energy(TRUE);

            if(write_stateFile){
              fprintf(stateFile,"\t<run_requested>%d</run_requested>\n",nruns);
              fprintf(stateFile,"\t<runs>\n");
            }
            for (j=0; j<nruns; j++) {

		Real eintra = 0.0;  // sum of intramolecular energy for the ligand plus that of the protein
		Real einter = 0.0; // intermolecular energy between the ligand and the protein
		struct tms tms_runStart, tms_runEnd;
		Clock  runStart, runEnd;
		FILE *tlogFile;
#ifdef _OPENMP
/* MPique TODO*/
		int thread_num=omp_get_thread_num();
		if(nruns>1) tlogFile = threadLogOpen( j );
		else tlogFile=logFile;
		if(tlogFile==NULL) stop("failed to create thread log file");
#ifdef DEBUG2
		fprintf(tlogFile, "run %2d nconf %2d GALS/GS on thread_num %d\n",
	          nconf+j+1, nconf, thread_num);
		fflush(tlogFile);
#endif
		/* MPique TODO add interrupt handler to remove all tlog files */
#else
		tlogFile=logFile;
#endif
		
		if(outlev>LOGBASIC) 
                (void) fprintf( tlogFile, "\n\tBEGINNING %s DOCKING %d of %d\n", 
		GlobalSearchMethod->longname(), j+1, nruns);

		/* set RNG seed using global run number, except for the first run in the job 
		 *   which continues on with seeds as possibly modified by ligand state
		 *   randomizations done during job setup but after "seed" DPF token,
		 *   (e.g. tran0/dihe0 or extended conformation search). The sole reason
		 *   for excepting the first run is for compatibility with existing tests.
		 *   The reason for the "getsd()" is so the actual first-run seeds will be
		 *   correctly reported so the run can be reproduced if necessary.
		 *   M Pique Oct 2013
		 */
		if(j==0&&nconf==0) getsd(&runseed[0][0], &runseed[0][1]);
		setsd(runseed[nconf+j][0], runseed[nconf+j][1]); 

                pr( tlogFile, "Run: %d Seed: %ld %ld [ Run %d of %d GA/GALS ]\n", nconf+j+1,
		 (long)runseed[nconf+j][0], (long)runseed[nconf+j][1],
		 j+1, nruns );

		if(outlev>=LOGRUNVV) {
			pr(tlogFile, "Date:\t");
			printdate( tlogFile, 2 );
			(void) fflush( tlogFile );
		}

                runStart = times( &tms_runStart );

		// MP TODO what does this do? should it be per-thread?
                if (B_reorient_random == TRUE) {
                    // create a new random orientation before docking
                    // reorient the ligand
                    reorient( tlogFile, true_ligand_atoms, atomstuff, crdpdb, charge, type,
                              parameterArray, randomQuat(), origin, ntor, tlist, vt, &ligand, debug, outlev );
                    // update the evaluate object
                    evaluate.update_crdpdb( crdpdb, vt );
                }

                //  Can get rid of the following line
                ((Genetic_Algorithm *)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor, outlev, logFile);

                // Reiterate output level...
		if(outlev>=LOGRUNV)
                pr(tlogFile, "Output level is set to %d.\n\n", outlev);

                // Start Lamarckian GA run -- Bound simulation
                sHist[nconf+j] = call_glss( GlobalSearchMethod, LocalSearchMethod,
                                          sInit,
                                          num_evals, pop_size,
                                          outlev, tlogFile,
                                          output_pop_stats, &ligand, &evaluate,
                                          B_RandomTran0, B_RandomQuat0, B_RandomDihe0,
                                          info, FN_pop_file, end_of_branch );
                // State of best individual at end of GA-LS run is returned.
                // Finished Lamarckian GA run
                
                runEnd = times( &tms_runEnd );
                if(outlev>=LOGRUNVV) {
			pr( tlogFile, "\nRun completed;  time taken for this run:\n");
			timesyshms( runEnd - runStart, &tms_runStart, &tms_runEnd, tlogFile);
			printdate( tlogFile, 1 );
			}
                if(outlev>=LOGRUNV) {

			pr(tlogFile, "Total number of Energy Evaluations: %u\n", evaluate.evals() );
			pr(tlogFile, "Total number of Generations:        %u\n", ((Genetic_Algorithm *)GlobalSearchMethod)->num_generations());
			}
		(void) fflush( tlogFile );

		if(outlev>LOGBASIC) {
                pr( tlogFile, "\n\n\tFINAL %s DOCKED STATE\n", GlobalSearchMethod->longname());
                pr( tlogFile,     "\t_______________________________________________\n\n\n" );
		}

                writePDBQT( nconf+j, runseed[nconf+j],  FN_ligand, dock_param_fn, lig_center,
                            sHist[nconf+j], ntor, &eintra, &einter, natom, atomstuff,
                            crd, peratomE,
                            charge, abs_charge, qsp_abs_charge,
                            ligand_is_inhibitor,
                            torsFreeEnergy,
                            vt, tlist, crdpdb, nonbondlist,
                            ad_energy_tables,
                            type, 
			    Nnb, Nnb_array, &group_energy, true_ligand_atoms,
			    B_calcIntElec,
                            map,
                            ignore_inter,
                            B_include_1_4_interactions, scale_1_4, parameterArray, unbound_internal_FE,
                            info, DOCKED, PDBQT_record, B_use_non_bond_cutoff, B_have_flexible_residues, ad4_unbound_model,
			    outlev, tlogFile);

                // See also "calculateEnergies.cc", switch(ad4_unbound_model)
                if (ad4_unbound_model == Unbound_Same_As_Bound) {
                    // Treat the unbound internal energy as the current internal energy
                    econf[nconf+j] = einter + torsFreeEnergy;
                }
                else econf[nconf+j] = eintra + einter + torsFreeEnergy - unbound_internal_FE;

                if(outlev>LOGBASIC) pr( tlogFile, UnderLine );
#ifdef _OPENMP
		if(nruns>1) threadLogClose(j);
#endif
            } // Next LGA run
	    nconf += nruns;
#ifdef _OPENMP
		if(nruns>1) {
			fprintf(stderr," concat %d log files...",nruns); fflush(stderr); /* MP debug */
			for(j=0;j<nruns;j++) {
				threadLogConcat(logFile, j);
				threadLogFree(j);
				}
			fprintf(stderr," done.\n"); fflush(stderr); /* MP debug */
		}
#endif

            if(write_stateFile){
               fprintf(stateFile,"\t</runs>\n");
               (void) fflush(stateFile);
            }
        break;

//______________________________________________________________________________

    case DPF_LS:
       // ls_run  |  do_local_only
	       nfields = sscanf(line, "%*s %d", &nruns); // optional way to specify nruns
            if ( nruns+nconf > MAX_RUNS ) {

               prStr( error_message, "%s:  ERROR: %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", programname, nruns+nconf, MAX_RUNS);
               stop( error_message );
           } 
	   if (LocalSearchMethod==NULL) {

               prStr(error_message, "%s:  ERROR:  You must use \"set_sw1\", \"set_psw1\" or \"set_pattern\" to create a Local Optimization object.\n", programname);
               stop(error_message);
            }
           exit_if_missing_elecmap_desolvmap_about("ls");

	    // set lig_center if not already set, use to center "crdpdb" ligand
	    center_ligand(crdorig, !B_found_about_keyword, natom, true_ligand_atoms,
	      tlist, ntor, crdpdb, lig_center, &sInit.Center, &ligand.S.Center,
	      outlev>=LOGBASIC, outlev, logFile);
	    // save centered crdpdb coords as crd (not sure is needed - MP 2012
	    for(int a=0;a<natom;a++) for(xyz=0;xyz<SPACE;xyz++)  
	      crd[a][xyz]=crdpdb[a][xyz];
	    // set 'tran0' vector to same as 'about' if not specified (2011-09)
	    if ( ! B_found_tran0_keyword ) {
		    ligand.S.T = sInit.T = sInit.Center;
		    pr( logFile, 
		    "Setting 'tran0' value to same as 'about' value: %.3f %.3f %.3f\n",
		    ligand.S.T.x, ligand.S.T.y, ligand.S.T.z);
		    B_found_tran0_keyword = TRUE;
	    }
           pr( logFile, "Number of Local Search (LS) only dockings = %d run%s\n", nruns, pl(nruns));
           if (ad4_unbound_model==Unbound_Default) ad4_unbound_model = Unbound_Same_As_Bound;
           pr(logFile, "Unbound model to be used is %s.\n", report_parameter_library());
/* MP experimental TODO  moved into per-thread code
           evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom,
                           info, map, peratomE,
                           nonbondlist,
                           ad_energy_tables,
                           Nnb, Nnb_array, &group_energy,
			   B_calcIntElec, B_isGaussTorCon,B_isTorConstrained,
                           B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, sInit, ligand,
                           ignore_inter,
                           B_include_1_4_interactions, scale_1_4, scale_eintermol,
                           unbound_internal_FE, B_use_non_bond_cutoff, B_have_flexible_residues, 
			   true_ligand_atoms, outlev, logFile);

            evaluate.compute_intermol_energy(TRUE);
*/

           if(write_stateFile){
             fprintf(stateFile,"\t<run_requested>%d</run_requested>\n",nruns);
             fprintf(stateFile,"\t<runs>\n");
           }

		{ // MP begin storage allocation for hacks...

	   // MP hack: create all necessary rho_ptrs, lb_rho_ptrs, tLocalSearchMethods
	   //  MP: Note this overrides the LocalSearchMethod defined in the DPF,
	   //  which might have been non-pseudo SW, or even "pattern"
	   //  I will fix this when the virtual functions are revised.
	   Real *trho_ptr[NUMG], *tlb_rho_ptr[NUMG];
	   Pseudo_Solis_Wets1 *tLocalSearchMethod[NUMG];  // should be clone of LocalSearchMethod object MP
	   Eval *tevaluate[NUMG];

	   for(int t=0;t<NUMG;t++) {
		trho_ptr[t] = new Real[7+sInit.ntor];
		tlb_rho_ptr[t] = new Real[7+sInit.ntor];
		tLocalSearchMethod[t] = 
                   new Pseudo_Solis_Wets1(7+sInit.ntor, max_its, max_succ, max_fail, 2.0, 0.5, trho_ptr[t], tlb_rho_ptr[t]);
 		tevaluate[t] = new Eval;
	   }
#pragma omp parallel for \
   shared(crdpdb,nconf) \
   private(j) \
   schedule(static)
           for (j=0; j<nruns; j++) {

		/* per-thread private locals: */
		EnergyComponent tperatomE[MAX_ATOMS];
		GroupEnergy tgroup_energy; // energy components of each of the five groups (intra-ligand, inter, and intra-receptor...)
		Real tcrd[MAX_ATOMS][SPACE];     // current coordinates according to State
		//Eval tevaluate;
		FILE *tlogFile;
		Real eintra = 0.0;  // sum of intramolecular energy for the ligand plus that of the protein
		Real einter = 0.0; // intermolecular energy between the ligand and the protein
		struct tms tms_runStart, tms_runEnd;
		Clock  runStart, runEnd;
		int tn; // thread number 0..NUMG
		//Pseudo_Solis_Wets1 *tLocalSearchMethod; // TODO clone of *LocalSearchMethod
		//MP not needed here   Real *rho_ptr ; // for PSW array of rho
		//MP not needed here   Real *lb_rho_ptr ; // for PSW array of lb_rho
		int d; // PSW variable index
		
		tn=omp_get_thread_num();

		if(nruns>1 && omp_get_max_threads()>1) tlogFile = threadLogOpen( j );
		else tlogFile=logFile;
		if(tlogFile==NULL) stop("failed to create thread log file");

	    // save centered crdpdb coords as crd (not sure is needed - MP 2012
	    for(int a=0;a<natom;a++) for(int xyz=0;xyz<SPACE;xyz++)  
	      tcrd[a][xyz]=crdpdb[a][xyz];

	   tevaluate[tn]->reset();
           tevaluate[tn]->setup( tcrd, charge, abs_charge, qsp_abs_charge, type, natom,
                           info, map, tperatomE,
                           nonbondlist,
                           ad_energy_tables,
                           Nnb, Nnb_array, &tgroup_energy,
			   B_calcIntElec, B_isGaussTorCon,B_isTorConstrained,
                           B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, sInit, ligand,
                           ignore_inter,
                           B_include_1_4_interactions, scale_1_4, scale_eintermol,
                           unbound_internal_FE, B_use_non_bond_cutoff, B_have_flexible_residues, 
			   true_ligand_atoms, outlev, tlogFile);

            /* this is the default:*/ tevaluate[tn]->compute_intermol_energy(TRUE);
		/* set RNG seed using global run number */
		if(nconf+j==0) getsd(&runseed[nconf][0], &runseed[nconf][1]);
		setsd(runseed[nconf+j][0], runseed[nconf+j][1]); 

	       if(outlev>LOGBASIC)
               (void) fprintf( tlogFile, "\tBEGINNING SOLIS & WETS LOCAL SEARCH DOCKING\n");
                pr( tlogFile, "Run: %d Seed: %ld %ld  [ Run %d of %d LS ]\n", nconf+j+1,
		 (long)runseed[nconf+j][0], (long)runseed[nconf+j][1],
		  j+1, nruns );

               pr(tlogFile, "Date:\t");
               printdate( tlogFile, 2 );
               (void) fflush( tlogFile );

               runStart = times( &tms_runStart );

 	       //tLocalSearchMethod = *LocalSearchMethod;  // copy and assign MP TODO 2014 not enough
		// start hacks here: MP TODO
      //  Allocate space for the variable rho's
      //  The lb_rho_ptr[] values are really const
      //rho_ptr = new Real[7+sInit.ntor];
      //lb_rho_ptr = new Real[7+sInit.ntor];

      //  Initialize the rho's corresponding to the translation
      //  0,1,2   x,y,z
      //  3,4,5,6 qx,qy,qz,qw
      //  7,...   tor1
//These scale values can be changed in the dpf, officially unsupported 4/2009
//Real psw_trans_scale = 1.0;//1 angstrom
//Real psw_rot_scale = 0.05;  //about 3 degrees, we think
//Real psw_tors_scale = 0.1; //about 6 degrees

      for (d=0; d<3; d++) {
         // d=0,1,2
         trho_ptr[tn][d] = rho * psw_trans_scale;// formerly trnStep0;
         tlb_rho_ptr[tn][d] = lb_rho * psw_trans_scale; //once trnStepFinal;
      }

      //  Initialize the rho's corresponding to the quaterion
      for (; d<7; d++) {
         // d=3,4,5,6
         trho_ptr[tn][d] = rho * psw_rot_scale;// formerly qtwStep0;
         tlb_rho_ptr[tn][d] = lb_rho * psw_rot_scale; //once qtwStepFinal;
      }

      //  Initialize the rho's corresponding to the torsions
      for (; d<7+sInit.ntor; d++) {
         // d=7,...
         trho_ptr[tn][d] = rho * psw_tors_scale;// formerly torStep0;
         tlb_rho_ptr[tn][d] = lb_rho * psw_tors_scale;//formerly torStepFinal;
      }

      //tLocalSearchMethod = new Pseudo_Solis_Wets1(7+sInit.ntor, max_its, max_succ, max_fail, 2.0, 0.5, rho_ptr, lb_rho_ptr);


               sHist[nconf+j] = call_ls(tLocalSearchMethod[tn], sInit, pop_size, &ligand,
		tevaluate[tn], outlev, tlogFile);

   // hacks here TODO MP
	// dumps core if you try...    delete tLocalSearchMethod;
	//delete [] rho_ptr;
	//delete [] lb_rho_ptr;
	
               pr(tlogFile, "There were %u Energy Evaluations.\n", tevaluate[tn]->evals());

	       if(outlev>=LOGRUNV) {
                  runEnd = times( &tms_runEnd );
                  pr( tlogFile, "Time taken for this Local Search (LS) run:\n");
                  timesyshms( runEnd - runStart, &tms_runStart, &tms_runEnd, tlogFile );
                  pr( tlogFile, "\n");
		  }

	       if(outlev>=LOGFORADT) {
                  pr( tlogFile, "\n\n\tFINAL LOCAL SEARCH DOCKED STATE\n" );
                  pr( tlogFile,     "\t_______________________________\n\n\n" );
		  }

               writePDBQT( nconf+j, runseed[nconf+j], FN_ligand, dock_param_fn, lig_center,
                           sHist[nconf+j], ntor, &eintra, &einter, natom, atomstuff,
                           tcrd, tperatomE,
                           charge, abs_charge, qsp_abs_charge,
                           ligand_is_inhibitor,
                           torsFreeEnergy,
                           vt, tlist, crdpdb, nonbondlist,
                           ad_energy_tables,
                           type, Nnb, Nnb_array, &tgroup_energy, true_ligand_atoms,
			   B_calcIntElec,
                           map,
                           ignore_inter,
                           B_include_1_4_interactions, scale_1_4, parameterArray, /*MP unbound_internal_FE*/0.,
                           info, DOCKED, PDBQT_record, B_use_non_bond_cutoff, B_have_flexible_residues, ad4_unbound_model,
			   outlev, tlogFile);

               // See also "calculateEnergies.cc", switch(ad4_unbound_model)
               if (ad4_unbound_model == Unbound_Same_As_Bound) {
                   // Treat the unbound internal energy as the current internal energy
                   econf[nconf+j] =  einter + torsFreeEnergy;
               }
               else econf[nconf+j] =  einter + eintra + torsFreeEnergy - unbound_internal_FE;

               if(outlev>=LOGRUNV) pr( tlogFile, UnderLine );
               (void) fflush( tlogFile );
		if(nruns>1 && omp_get_max_threads()>1) threadLogClose( j );

           } // Next run - also close of 'parallel for' region
	   // MP hack:
	   for(int t=0;t<NUMG;t++) {
		if(tevaluate[t]!=NULL) delete tevaluate[t];
		if(tLocalSearchMethod[t]!=NULL) delete tLocalSearchMethod[t];
		// done by PSWLocalSearchMethod destructor:  delete [] trho_ptr[t];
		// done by PSWLocalSearchMethod destructor:  delete [] tlb_rho_ptr[t];
	   }
       	} // MP end storage allocation for hacks...
	   nconf += nruns;
		if(nruns>1 && omp_get_max_threads()>1) {
			fprintf(stderr," concat log files..."); fflush(stderr); /* MP debug */
			for(int j=0;j<nruns;j++) {
				threadLogConcat(logFile, j);
				threadLogFree(j);
				}
			fprintf(stderr," done.\n"); fflush(stderr); /* MP debug */
		}
           if(write_stateFile){
             fprintf(stateFile,"\t</runs>\n");
             (void) fflush(stateFile);
           }
       break;

//______________________________________________________________________________

    case GA_pop_size:
       get1arg(line, "%*s %u", &pop_size, "GA_POP_SIZE");
       pr(logFile, "A population of %u individuals will be used\n", pop_size);
       break;

//______________________________________________________________________________

    case GA_num_generations:
       get1arg(line, "%*s %u", &num_generations, "GA_NUM_GENERATIONS");
       pr(logFile, "The GA will run for at most %u generations.\n", num_generations);
       break;

//______________________________________________________________________________

    case GA_num_evals:
       get1arg(line, "%*s %u", &num_evals, "GA_NUM_EVALS");
       pr(logFile, "There will be at most %u function evaluations used.\n", num_evals);
       break;

//______________________________________________________________________________

    case GA_window_size:
       get1arg(line, "%*s %d", &window_size, "GA_WINDOW_SIZE");
       pr(logFile, "The GA's selection window is %d generations.\n", window_size);
       break;

//______________________________________________________________________________

    case GA_low:
       get1arg(line, "%*s %d", &low, "GA_LOW");
       pr(logFile, "Setting GA low to %d.\n", low);
       break;

//______________________________________________________________________________

    case GA_high:
       get1arg(line, "%*s %d", &high, "GA_HIGH");
       pr(logFile, "Setting GA high to %d.\n", high);
       break;

//______________________________________________________________________________

    case GA_elitism:
       get1arg(line, "%*s %d", &elitism, "GA_ELITISM");
       pr(logFile, "The %d best will be preserved each GA generation.\n", elitism);
       break;

//______________________________________________________________________________

    case GA_mutation_rate:
       get1arg(line, "%*s " FDFMT, &m_rate, "GA_MUTATION_RATE");
       pr(logFile, "The mutation rate is %f.\n", m_rate);
      // if m_rate is out of range, make_table will fail 
      if (m_rate < 0 || m_rate > 1) {
	prStr(error_message, "mutation rate must be within range 0 to 1 (inclusive).\n");
	stop(error_message);
	}
       break;

//______________________________________________________________________________

    case GA_crossover_rate:
       get1arg(line, "%*s " FDFMT, &c_rate, "GA_CROSSOVER_RATE");
       pr(logFile, "The crossover rate is %f.\n", c_rate);
       break;

//______________________________________________________________________________

    case GA_Cauchy_alpha:
       get1arg(line, "%*s " FDFMT, &alpha, "GA_CAUCHY_ALPHA");
       pr(logFile, "The alpha parameter (for the Cauchy distribution) is being set to %f.\n",
          alpha);
       break;

//______________________________________________________________________________

    case GA_Cauchy_beta:
       get1arg(line, "%*s " FDFMT, &beta, "GA_CAUCHY_BETA");
       pr(logFile, "The beta parameter (for the Cauchy distribution) is being set to %f.\n",
          beta);
       break;

//______________________________________________________________________________

    case SW_max_its:
       get1arg(line, "%*s %u", &max_its, "SW_MAX_ITS");
       pr(logFile, "Solis & Wets algorithms will perform at most %u iterations.\n", max_its);
       break;

//______________________________________________________________________________

    case SW_max_succ:
       get1arg(line, "%*s %u", &max_succ, "SW_MAX_SUCC");
       pr(logFile, "Solis & Wets algorithms expand rho every %u in a row successes.\n", max_succ);
      break;

//______________________________________________________________________________

    case SW_max_fail:
       get1arg(line, "%*s %u", &max_fail, "SW_MAX_FAIL");
       pr(logFile, "Solis & Wets algorithms contract rho every %u in a row failures.\n", max_fail);
      break;

//______________________________________________________________________________

    case SW_rho:
       get1arg(line, "%*s " FDFMT, &rho, "SW_RHO");
       pr(logFile, "rho is set to %f.\n", rho);
      break;


//______________________________________________________________________________

    case SW_lb_rho:
        get1arg(line, "%*s " FDFMT, &lb_rho, "SW_LB_RHO");
        pr(logFile, "rho will never get smaller than %f.\n", lb_rho);
        break;

//______________________________________________________________________________

    case PSW_TRANS_SCALE:
        get1arg(line, "%*s " FDFMT, &psw_trans_scale, "PSW_TRANS_SCALE");
        pr(logFile, "psw_trans_scale is set to %f.\n", psw_trans_scale);
        break;

//______________________________________________________________________________

    case PSW_ROT_SCALE:
        get1arg(line, "%*s " FDFMT, &psw_rot_scale, "PSW_ROT_SCALE");
        pr(logFile, "psw_rot_scale is set to %f.\n", psw_rot_scale);
        break;

//______________________________________________________________________________

    case PSW_TORS_SCALE:
        get1arg(line, "%*s " FDFMT, &psw_tors_scale, "PSW_TORS_SCALE");
        pr(logFile, "psw_tors_scale is set to %f.\n", psw_tors_scale);
        break;

//______________________________________________________________________________


    case LS_search_freq:
        get1arg(line, "%*s " FDFMT, &localsearch_freq, "LS_SEARCH_FREQ");
        pr(logFile, "Local search will be performed with frequency %f.\n", localsearch_freq);
        break;

//______________________________________________________________________________
   
   case PSO_C1:
        get1arg(line, "%*s %lf", &pso_options.c1, "PSO_C1");
        pr(logFile, "PSO will be performed with the First Coefficient  (C1) %lf.\n", pso_options.c1);
        break;

//______________________________________________________________________________

   case PSO_C2:
        get1arg(line, "%*s %lf", &pso_options.c2, "PSO_C2");
        pr(logFile, "PSO will be performed with the Second Coefficient (C2) %lf.\n", pso_options.c2);
        break;

//______________________________________________________________________________

   case PSO_K:
        get1arg(line, "%*s %d", &pso_options.pso_K, "PSO_K");
	// MP TODO this should be allowed to be equal to pop_size
	if(pso_options.pso_K>PSO_K_MAX) stop("PSO_K_max too big ");
        pr(logFile, "Max number of PSO particles informing a given one = %d.\n", pso_options.pso_K);
        break;

//______________________________________________________________________________
//
//   case PSO_swarm_moves:
//        get1arg(line, "%*s %d", &eval_max, "PSO_SWARM_MOVES");
//        pr(logFile, "There will be %d swarm Moves.\n", eval_max);
//        break;
//
//______________________________________________________________________________

//   case PSO_swarm_size_factor:
//        get1arg(line, "%*s %d", &S_factor, "PSO_SS_FACTOR");
//        pr(logFile, "There will be %d Swarm Size Factor.\n", S_factor);
//        break;
//______________________________________________________________________________

//   case PSO_n_exec:
//        get1arg(line, "%*s %d", &n_exec_max, "PSO_N_EXEC");
//        pr(logFile, "Number of requested PSO runs = %d.\n", n_exec_max);
//        break;

//______________________________________________________________________________
        
 	case PSO_W_START:
	     get1arg( line, "%*s %lf", &pso_options.pso_w_start, "PSO_W_START");
        pr(logFile, "PSO starting velocity weight  = %f.\n", pso_options.pso_w_start);
		break;
//______________________________________________________________________________
 
	case PSO_W_END:
	    get1arg( line, "%*s %lf", &pso_options.pso_w_end, "PSO_W_END");
        pr(logFile, "PSO ending velocity weight  = %f.\n", pso_options.pso_w_end);
		break;
//______________________________________________________________________________
 
	case PSO_NEIGHBORS_DYNAMIC:
	    get1arg( line, "%*s %d", &pso_options.pso_neighbors_dynamic,
	     "PSO_NEIGHBORS_DYNAMIC");
        pr(logFile, "PSO neighbors dynamic  = %d.\n", mkbool(pso_options.pso_neighbors_dynamic));
		break;
//______________________________________________________________________________
 
	case PSO_NEIGHBORS_SYMMETRIC:
	    get1arg( line, "%*s %d", &pso_options.pso_neighbors_symmetric,
	     "PSO_NEIGHBORS_SYMMETRIC");
        pr(logFile, "PSO neighbors dynamic  = %d.\n", mkbool(pso_options.pso_neighbors_dynamic));
		break;
//______________________________________________________________________________
 
	case PSO_RANDOM_BY_DIMENSION:
	    get1arg( line, "%*s %d", &pso_options.pso_random_by_dimension, "PSO_RANDOM_BY_DIMENSION");
        pr(logFile, "PSO random by dimension  = %d.\n", mkbool(pso_options.pso_random_by_dimension));
		break;
//______________________________________________________________________________
 
	case PSO_ADAPTIVE_VELOCITY:
	    get1arg( line, "%*s %d", &pso_options.pso_adaptive_velocity,
	     "PSO_ADAPTIVE_VELOCITY");
        pr(logFile, "PSO adaptive velocity  = %d.\n", mkbool(pso_options.pso_adaptive_velocity));
		break;
//______________________________________________________________________________

	case PSO_REGENERATE_AT_LIMIT:
	    get1arg( line, "%*s %d", &pso_options.pso_regenerate_at_limit,
	     "REGENERATE_AT_LIMIT");
        pr(logFile, "PSO regenerate_at_limit = %d.\n", mkbool(pso_options.pso_regenerate_at_limit));
		break;
//______________________________________________________________________________
 
	case PSO_STAGE2CONSTRICTION:
	    get1arg( line, "%*s %d", &pso_options.pso_stage2constriction,
	     "PSO_STAGE2CONSTRICTION");
        pr(logFile, "PSO stage2constriction  = %d.\n", mkbool(pso_options.pso_stage2constriction));
		break;
//______________________________________________________________________________
 
	case PSO_INTERPOLATE_AS_SCALARS:
	    get1arg( line, "%*s %d", &pso_options.pso_interpolate_as_scalars,
	     "PSO_INTERPOLATE_AS_SCALARS");
        pr(logFile, "PSO interpolate as scalars = %d.\n", mkbool(pso_options.pso_interpolate_as_scalars));
		break;
// 	
//	case PSO_OUTPUT_GENS:
//	    (void) sscanf( line, "%*s %d", &pso_output_gens);	    
//		break;
//_______________________________________________________________________________________
// parameters for pso vmax. vmin = -vmax
//_______________________________________________________________________________________	    
	case PSO_TVMAX:
	    get1arg( line, "%*s %f", &pso_tvmax, "PSO_TVMAX");
		break;
//_______________________________________________________________________________________	    
	case PSO_QVMAX:
	    get1arg( line, "%*s %f", &pso_qvmax, "PSO_QVMAX");
		break;
//_______________________________________________________________________________________	    
	case PSO_RVMAX:
	    get1arg( line, "%*s %f", &pso_rvmax, "PSO_RVMAX");
		pso_rvmax = DegreesToRadians(pso_rvmax);
		break;	    
//_______________________________________________________________________________________	    
// 201105050 START
//______________________________________________________
////////////////////////////////////////////////////////////////////////////////
 // Entry point to call constriction CPSO
 ////////////////////////////////////////////////////////////////////////////////
    case DPF_PARSWARMOPT:
        int D; //search space dimension 
    	get1arg( line, "%*s %d", &nruns, "PARSWARMOPT" );
        
      //set GlobalSearchMethod to ParticleSwarmGS
	  if (GlobalSearchMethod != NULL) {
	 	pr(logFile, "Deleting the previous settings for the PSO.\n");	          
	    delete GlobalSearchMethod;
	    GlobalSearchMethod = NULL;
	  }

	 	pr( logFile, "\nTotal number of torsions in system = %d \n", sInit.ntor);
	 	D = 7 + ntor; //Dimension D of degree of freedom for a ligand
	 	pr(logFile, "\nTotal number of dimension is equal to the number of Degrees of Freedom = %d\n", D);	 
	    pr( logFile, "Number of requested PSO dockings = %d run%s\n", nruns, pl(nruns));
	 		 
	 	//No. of Particles, size of the Swarm 
	 	S = pop_size; //S - Swarm Size
	 	pr( logFile, "\nTotal number of particles in swarm = %d\n", S);	          	
	 	fflush(logFile);  

        
	  
	  if(outlev>LOGBASIC) {
		pr(logFile, "PSO tvmax = %.3f\n", pso_tvmax);
		pr(logFile, "PSO qvmax = %.3f\n", pso_qvmax);
		pr(logFile, "PSO rvmax = %.3f\n", pso_rvmax);
		}

	  float pso_vmax [PSO_D_MAX]; // TODO = float [7 + sInit.ntor];
	  float pso_vmin [PSO_D_MAX]; // TODO = float [7 + sInit.ntor];
	  //pso_vmax = float [7 + sInit.ntor];
	  //pso_vmin = float [7 + sInit.ntor];
	  //translation, quaternion, torsion
	  	//translation x, y, z
	  	for(j= 0; j < 3; j++) {
			pso_vmax[j] = pso_tvmax;
		  	pso_vmin[j] = -pso_tvmax;
		}
		// quaternion components  MP this isnt meaningful should be just a scalar angle
		for(j=3; j < 7; j++) {		
			pso_vmax[j] = pso_qvmax;	
			pso_vmin[j] = -pso_qvmax;			
		}						 			 	
	  // torsion part	  
	  for( j=7; j < 7+ sInit.ntor; j++) {
	  	pso_vmax[j] = pso_rvmax;	
		pso_vmin[j] = -pso_rvmax;		
        }
      
        if (outlev > 2) {
            pr( logFile, "\ncalling PSO initialiseDimension: info->lo=%f %f %f , info->hi=%f %f %f,  \n", info->lo[0], info->lo[1], info->lo[2], info->hi[0], info->hi[1], info->hi[2]);
	    };
        initialiseDimension(info, pso_xmin, pso_xmax, D);
        if (outlev > 2) {
            pr( logFile, "\nAFTER PSO initDim trans: xmin=%f %f %f , xmax=%f %f %f,  \n", pso_xmin[0], pso_xmin[1], pso_xmin[2], pso_xmax[0], pso_xmax[1], pso_xmax[2]);
            pr( logFile, "\nAFTER PSO initDim quat: xmin=%f %f %f %f, xmax=%f %f %f %f,  \n", pso_xmin[3], pso_xmin[4], pso_xmin[5], pso_xmin[6], pso_xmax[3], pso_xmax[4], pso_xmax[5],pso_xmax[6]);
        };
	  

	  GlobalSearchMethod = new ParticleSwarmGS(
                        pso_vmax,
                        pso_vmin,
                        pso_xmax,
                        pso_xmin,
						pso_options,
	  					LocalSearchMethod, 
						num_evals,
	  					num_generations, 
	  					output_pop_stats);  
        ((ParticleSwarmGS*)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor, outlev, logFile);
	  					   
      pr(logFile, "GlobalSearchMethod is set to PSO.\n\n");     
   	 

	    // set lig_center if not already set, use to center "crdpdb" ligand
	    center_ligand(crdorig, !B_found_about_keyword, natom, true_ligand_atoms,
	      tlist, ntor, crdpdb, lig_center, &sInit.Center, &ligand.S.Center,
	      outlev>=LOGBASIC, outlev, logFile);
	    // save centered crdpdb coords as crd (not sure is needed - MP 2012
	    for(int a=0;a<natom;a++) for(xyz=0;xyz<SPACE;xyz++)  
	      crd[a][xyz]=crdpdb[a][xyz];

	   evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom,
                        info, map, peratomE, nonbondlist, ad_energy_tables, Nnb,
			Nnb_array, &group_energy,
                        B_calcIntElec, B_isGaussTorCon, B_isTorConstrained,
                        B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, sInit, ligand,
                        ignore_inter, B_include_1_4_interactions, scale_1_4, scale_eintermol,
                        unbound_internal_FE, B_use_non_bond_cutoff, B_have_flexible_residues,
			true_ligand_atoms, outlev, logFile);
	 		 		 	
        evaluate.compute_intermol_energy(TRUE);
   
	   //BEGINNING PARTICLE SWARM OPTIMIZATION run
	   for (j = 0; j < nruns; j++)
	   {	 
		Real eintra = 0.0;  // sum of intramolecular energy for the ligand plus that of the protein
		Real einter = 0.0; // intermolecular energy between the ligand and the protein
		struct tms tms_runStart, tms_runEnd;
		Clock  runStart, runEnd;
		//(void) fprintf( logFile, "\n\tBEGINNING PARTICLE SWARM OPTIMIZATION (PSO) \n");
	    if(outlev>LOGBASIC)
            (void) fprintf( logFile, "\n\tBEGINNING %s DOCKING\n", GlobalSearchMethod->longname());

		/* set RNG seed using global run number */
		if(nconf==0&&j==0) getsd(&runseed[nconf][0], &runseed[nconf][1]);
		else setsd(runseed[nconf+j][0], runseed[nconf+j][1]); 

                pr( logFile, "Run: %d Seed: %ld %ld [ Run %d of %d %s ]\n", nconf+j+1,
		 (long)runseed[nconf+j][0], (long)runseed[nconf+j][1],
		 j+1, nruns, GlobalSearchMethod->shortname() );
	 		pr(logFile, "Date:\t");
	        printdate(logFile, 2 );
	 		(void)fflush(logFile);
		 										
	 		runStart = times(&tms_runStart);
	 		//pr( logFile, "\nTotal number of torsions in system = %d \n", sInit.ntor);
	 		//Start Particle Swarm Optimization Run	               	 		
                sHist[nconf+j] = call_glss( GlobalSearchMethod, LocalSearchMethod,
                                          sInit,
                                          num_evals, pop_size,
                                          outlev, logFile,
                                          output_pop_stats, &ligand, &evaluate,
                                          B_RandomTran0, B_RandomQuat0, B_RandomDihe0,
                                          info, FN_pop_file, end_of_branch );
	 		//Finished Particle Swarm Optimization Run	 	 		
	 		runEnd = times(&tms_runEnd);
            pr(logFile, "Time taken for this PSO run:\n");
            timesyshms(runEnd-runStart, &tms_runStart, &tms_runEnd, logFile);
            pr(logFile, "\n");
            (void) fflush(logFile);
	 				
	 				        
	        pr(logFile, "Total number of Energy Evaluations: %u\n", evaluate.evals());	        	        
            pr(logFile, "Total number of Generations:        %u\n", ((ParticleSwarmGS*)GlobalSearchMethod)->num_generations()); // TSRI 20101101 added by M Pique
	 							 		
	 		pr( logFile, "\n\n\tFINAL PSO DOCKED STATE\n" );
	 		pr( logFile, "\t____________________________________________________________\n\n\n" );
	 		
             writePDBQT( nconf+j, runseed[nconf],  FN_ligand, dock_param_fn, lig_center,
                       sHist[nconf], ntor, &eintra, &einter, natom, atomstuff,
                       crd, peratomE, charge, 
                       abs_charge, qsp_abs_charge,
                       ligand_is_inhibitor,
                       torsFreeEnergy,
                       vt, tlist, crdpdb, nonbondlist,
                       ad_energy_tables,
                       type, 
		       Nnb, Nnb_array, &group_energy, true_ligand_atoms,
		       B_calcIntElec, map,
                       ignore_inter, B_include_1_4_interactions, 
                       scale_1_4, parameterArray, unbound_internal_FE,
                       info, DOCKED, PDBQT_record, B_use_non_bond_cutoff, //info
                       B_have_flexible_residues, ad4_unbound_model,
		       outlev, logFile);
       	

            //econf[nconf+j] = eintra + einter; // changed to next line M Pique June 2013
                  econf[nconf+j] = eintra + einter + torsFreeEnergy - unbound_internal_FE;
	 		pr( logFile, UnderLine );	
	 		 			 			 				
	 	} // next PSO run j
	    nconf += nruns;
	 	
	 	if(write_stateFile){
           fprintf(stateFile,"\t</runs>\n");
           (void) fflush(stateFile);
        }
        (void) fflush(logFile);          	 	 
	    break;
	    // end  DPF_CPSO_RUN
//________________________________________________________________________
// 20110505 END 
//________________________________________________________________________

    case DPF_ANALYSIS:
        /*
        ** analysis
        */
        /* _____________________________________________________________________
        **
        ** Perform Cluster analysis on results of docking,
	** across all "nconf" runs so far, which might be a series of local, global,
	** simanneal, or hybrid runs
        ** _____________________________________________________________________
        */
            analysis( Nnb, Nnb_array, &group_energy, true_ligand_atoms,
	              atomstuff, charge, abs_charge, qsp_abs_charge, B_calcIntElec, clus_rms_tol,
                      crdpdb, ad_energy_tables, map, econf, nconf,
                      natom, nonbondlist, nconf, ntor, sHist, FN_ligand,
                      lig_center, B_symmetry_flag, B_unique_pair_flag, tlist, type, vt, FN_rms_ref_crds,
                      torsFreeEnergy, B_write_all_clusmem, ligand_is_inhibitor,
                      ignore_inter, B_include_1_4_interactions, scale_1_4,
                      unbound_internal_FE,
                      info, B_use_non_bond_cutoff, B_have_flexible_residues,
                      B_rms_atoms_ligand_only, ad4_unbound_model, 
                      B_rms_heavy_atoms_only, h_index, outlev, logFile);

        break;

//______________________________________________________________________________

    case DPF_TORSDOF:
        /*
        ** torsdof %d %f
        */
        nfields = sscanf( line, "%*s %d " FDFMT, &ntorsdof, &torsdoffac );
        if (nfields == 2) {
            pr( logFile, "WARNING:  The torsional DOF coefficient is now read in from the parameter file; the value specified here (%.4lf) will be ignored.\n\n", (double)torsdoffac);
	// TODO should this be fatal error?  M Pique 2010
        }
        pr( logFile, "Number of torsional degrees of freedom = %d\n", ntorsdof);
        pr( logFile, "Free energy coefficient for torsional degrees of freedom = %.4f", AD4.coeff_tors);
        if (parameter_library_found) {
            pr( logFile, " as specified in parameter library \"%s\".\n\n", FN_parameter_library );
        } else {
            pr( logFile, ", the factory default value.\n\n");
        }

        torsFreeEnergy = (Real)ntorsdof * AD4.coeff_tors;

        pr( logFile, "Estimated loss of torsional free energy upon binding = %+.4f kcal/mol\n\n", torsFreeEnergy);
        break;

//______________________________________________________________________________

    case DPF_INVESTIGATE:
        /*
        ** Bin energies by RMSD from reference structure
        **
        ** investigate 100000 1000000 100
        */
        nfields = sscanf( line, "%*s %d %d %d", &OutputEveryNTests, &maxTests, &NumLocalTests );
	if(nfields!=3) stop("syntax error in INVESTIGATE or BIN_ENERGIES_BY_RMSD line");
        (void) fprintf( logFile, "OutputEveryNTests= %d\n", OutputEveryNTests);
        (void) fprintf( logFile, "maxTests= %d\n", maxTests );
        (void) fprintf( logFile, "NumLocalTests= %d\n\n", NumLocalTests );
	// M Pique TODO this probably should not use B_unique_pair_flag 2010
        (void) investigate( Nnb, Nnb_array, &group_energy,
	                    charge, abs_charge, qsp_abs_charge, B_calcIntElec,
                            crd, crdpdb, ad_energy_tables,
                            maxTests,
                            map, natom, nonbondlist, ntor,
                            tlist, type, vt, B_isGaussTorCon, US_torProfile,
                            B_isTorConstrained, B_ShowTorE, US_TorE,
                            F_TorConRange, N_con, B_symmetry_flag, B_unique_pair_flag, FN_rms_ref_crds,
                            OutputEveryNTests, NumLocalTests, trnStep0, torStep0,
                            ignore_inter,
                            B_include_1_4_interactions, scale_1_4, scale_eintermol,
                            unbound_internal_FE,
                            info, B_use_non_bond_cutoff, B_have_flexible_residues, 
                            B_rms_heavy_atoms_only, h_index,
			    true_ligand_atoms, outlev, logFile);

        break;

//______________________________________________________________________________

    case DPF_LIG_NOT_INHIB:
        /*
        ** ligand_is_not_inhibitor
        */
        ligand_is_inhibitor = 0;
        pr( logFile, "\nThis ligand is not an inhibitor, so dissociation constants (Kd) will be calculated, not inhibition constants (Ki).\n\n" );
        break;

/*____________________________________________________________________________*/

    case DPF_UNBOUND_MODEL:
        /*
        **  unbound_model { extended [energy <FLOAT>]| compact | bound }
        **    extended is alias for "compute_unbound_extended" token
        */
        char unbound_model_type[LINE_LEN];
        get1arg( line, "%*s %s", unbound_model_type, "UNBOUND_MODEL_TYPE" );

        if (streq( unbound_model_type, "bound")
        || streq( unbound_model_type, "same_as_bound")
        || streq( unbound_model_type, "unbound_same_as_bound")) {
            if (ad4_unbound_model != Unbound_Same_As_Bound)  // default for Autodock 4.1
                setup_parameter_library(logFile, outlev, "Unbound_Same_As_Bound", Unbound_Same_As_Bound, &AD4);
            ad4_unbound_model = Unbound_Same_As_Bound;
        } else if (streq( unbound_model_type, "extended")) {
            if (ad4_unbound_model != Unbound_Default) { //illegal to set extended after other
                pr( logFile, "%s:  ERROR:  Setting unbound model type twice: \"%s\" .\n",
                    programname, unbound_model_type );
                stop("");
            }
            if ( (1== sscanf( line, "%*s extended energy " FDFMT, &unbound_internal_FE ))){
                ad4_unbound_model = Extended;
                setup_parameter_library(logFile, outlev, "unbound_extended", ad4_unbound_model, &AD4);
            }
            else goto process_DPF_COMPUTE_UNBOUND_EXTENDED; // case DPF_COMPUTE_UNBOUND_EXTENDED below
        } else if (streq( unbound_model_type, "compact")) {
            ad4_unbound_model = Compact;
        } else {
            // note that "User" is not acceptable in dpf file
            pr( logFile, "%s:  ERROR:  Unrecognized unbound model type \"%s\" .\n",
                    programname, unbound_model_type );
            stop("");
        }
        break;

/*____________________________________________________________________________*/

    case DPF_UNBOUND:
        /*
         * unbound FLOAT
         * unbound energy FLOAT
         */
        if (ad4_unbound_model != Unbound_Default && ad4_unbound_model!= User) { //illegal to set user after other
            pr( logFile, "%s:  ERROR:  Setting unbound model type twice!\n",
                programname );
            stop("");
        }
        if ((1!= sscanf( line, "%*s " FDFMT, &unbound_internal_FE ))
        && (1!= sscanf( line, "%*s energy" FDFMT, &unbound_internal_FE ))){
            pr( logFile, "%s:  ERROR:  Non-numeric unbound model energy \"%s\" .\n",
                    programname, line);
            stop("Non-numeric unbound model energy");
        }
        pr(logFile, "The internal energy of the unbound state was set to %+.3lf kcal/mol\n", unbound_internal_FE);
        ad4_unbound_model = User;
        pr(logFile, "The unbound ligand energy model was set to User\n\n");
        break;

/*____________________________________________________________________________*/
    case DPF_COMPUTE_UNBOUND_EXTENDED:
        /*
         *  compute_unbound_extended
         */
    process_DPF_COMPUTE_UNBOUND_EXTENDED:
        if (ntor > 0) {
           (void) sprintf( message, "%s: WARNING: Using autodock4.0 unbound extended model in autodock4.2!\n", programname );
            print_2x( logFile, stderr, message );
            if (ad4_unbound_model != Unbound_Default) { //illegal to set extended after other
                pr( logFile, "%s:  ERROR:  Setting unbound model type twice!\n",
                    programname );
                stop("");
            }
            ad4_unbound_model = Extended;
            setup_parameter_library(logFile, outlev, "unbound_extended", ad4_unbound_model, &AD4);

            pr(logFile, "Computing the energy of the unbound state of the ligand,\ngiven the torsion tree defined in the ligand file.\n\n");

            // The initial goal is to obtain an extended conformation of the ligand.

            // Step 0 // {
            //
            // Set termination criteria for unbound calculations
            // -------------------------------------------------
            //
            // Set the maximum number of energy evaluations for finding the unbound conformation
            // if num_evals is less than this, then a shorter unbound docking will be performed
            max_evals_unbound = 1000000; // 1 million
            num_evals_unbound = num_evals > max_evals_unbound ?  max_evals_unbound :  num_evals;
            // end of Step 0 // }

            // Step 1 // {
            //
            // Run a hybrid global-local search using the unbound energy tables (set to be repulsive-only)
            // -------------------------------------------------------------------------------------------
            //
            //  *  Turn off the use of the non-bond cutoff
            //  *  Turn off internal electrostatics
            //  *  Turn off intermolecular energy calculations
            // TODO Need not translate or rotate the ligand in unbound searches
            //
            /*
             *  Genetic Algorithm-Local search,  a.k.a.
             *  Lamarckian Genetic Algorithm
             */
            if ((GlobalSearchMethod==NULL)||(LocalSearchMethod==NULL)) {
                prStr(error_message, "%s:  ERROR:  You must use \"set_ga\" to allocate both Global Optimization object AND Local Optimization object.\n", programname);
                stop(error_message);
            }
            exit_if_missing_elecmap_desolvmap_about("compute_unbound_extended");

            //
            // Do not use a non-bond cutoff, this helps to produce the "most" extended conformation
            // especially with long inhibitors
            B_use_non_bond_cutoff = FALSE;
            //
            // Save the current value of B_calcIntElec, so we can restore it later.
            B_calcIntElec_saved = B_calcIntElec;
            //
            // Set the calculation of internal electrostatics to FALSE:
            // B_calcIntElec = FALSE;
            //
            // Assume the unbound state of the receptor is the same as the input coordinates from the
            // flexible residues file.  This means we must not change the rotatable bonds in the
            // flexible residues of the receptor during the unbound extended search.
            // We can turn off rotation of the flexres by setting ntor to ntor_ligand.
            // Save the current value of "ntor" in the "sInit" state variable, set it to number of torsions
            // in the ligand for the unbound extended search, then restore it afterwards.
            saved_sInit_ntor = sInit.ntor;
            sInit.ntor = ntor_ligand;
            //
            // Use the repulsive unbound energy tables, "unbound_energy_tables",
            // to drive the molecule into an extended conformation
            evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom,
                            info, map, peratomE, nonbondlist, unbound_energy_tables, Nnb,
			    Nnb_array, &group_energy,
                            B_calcIntElec, B_isGaussTorCon, B_isTorConstrained,
                            B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, sInit, ligand,
                            ignore_inter,
                            B_include_1_4_interactions, scale_1_4, scale_eintermol,
                            unbound_internal_FE, B_use_non_bond_cutoff, B_have_flexible_residues,
			    true_ligand_atoms, outlev, logFile);
                            //parameterArray, unbound_internal_FE, B_use_non_bond_cutoff, B_have_flexible_residues);
            //
            // Turn off computing the intermolecular energy, we will only consider the intramolecular energy
            // to determine the unbound state of the flexible molecule:
            evaluate.compute_intermol_energy(FALSE);
            //
            (void) fprintf( logFile, "\n\tBEGINNING COMPUTATION OF UNBOUND EXTENDED STATE USING LGA\n");
            (void) fprintf( logFile,     "\t_________________________________________________________\n\n\n");
            //
            pr(logFile, "Date:\t");
            printdate( logFile, 2 );
            (void) fflush( logFile );
            //  Can get rid of the following line
            ((Genetic_Algorithm *)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor, outlev, logFile);
            //
            // Start Lamarckian GA run searching only torsions -- Unbound simulation
            // sUnbound_ext = call_glss_tors( GlobalSearchMethod, LocalSearchMethod,
            sUnbound_ext = call_glss( GlobalSearchMethod, LocalSearchMethod,
                                      sInit,
                                      num_evals_unbound, pop_size,
                                      outlev, logFile,
                                      output_pop_stats, &ligand, &evaluate,
                                      // B_RandomDihe0, // use this line with call_glss_tors()
                                      B_RandomTran0, B_RandomQuat0, B_RandomDihe0,
                                      info, FN_pop_file, end_of_branch );
            // State of best individual at end of GA-LS run, sUnbound_ext, is returned.
            // Finished Lamarckian GA run
            pr( logFile, "\nFinished Lamarckian Genetic Algorithm (LGA)\n");
            printdate( logFile, 1 );
            (void) fflush( logFile );
            pr(logFile, "\nTotal number of Energy Evaluations: %u\n", evaluate.evals() );
            pr(logFile, "Total number of Generations:        %u\n", ((Genetic_Algorithm *)GlobalSearchMethod)->num_generations());
            // end of Step 1 // }

            // Step 2 // {
            //
            // Do a short local search using the standard internal energy tables
            // -----------------------------------------------------------------
            //
            // turn on internal electrostatics
            // but keep intermolecular energy calculations off
            //
            // Turn on calculation of internal electrostatics:
            //// B_calcIntElec = TRUE;
            //
            // Use the standard AutoDock energy tables to compute the internal energy
            // Use this value to set unbound_internal_FE
            //// evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom, info, map,
                            //// peratomE, nonbondlist, ad_energy_tables, Nnb,
                            //// B_calcIntElec, B_isGaussTorCon, B_isTorConstrained,
                            //// B_ShowTorE, US_TorE, US_torProfile, vt, tlist, sInit, ligand,
                            //// ignore_inter,
                            //// B_include_1_4_interactions, scale_1_4,
                            //// unbound_internal_FE, B_use_non_bond_cutoff, B_have_flexible_residues,
            //
            // --- Start Local Search ---
            //// pr( logFile, "\nPerforming local search using standard AutoDock scoring function\n" );
            //// pr( logFile, "\nUsing UnboundLocalSearchMethod = new Solis_Wets1(7+sInit.ntor, 300, 4, 4, 1., 0.01, 2., 0.5, 1.);\n\n" );
            // Create a local search object
            // * Use an initial rho value of 0.1 (default is set in DPF by "sw_rho 1.0")
            //   to ensure smaller, 'more local' steps.
            // * Use a search frequency of 1.0 (default is set in DPF by "ls_search_freq 0.06")
            //// unsigned int ls_pop_size = 150;
            // max_its = 300
            // max_succ = 4
            // max_fail = 4
            // rho = 1.
            // lb_rho = 0.01
            // expansion = 2.
            // contraction = 0.5
            // search_freq = 1.
            //// UnboundLocalSearchMethod = new Solis_Wets1(7+sInit.ntor, 300, 4, 4, 1., 0.01, 2., 0.5, 1.);
            // Perform a local search, using the standard AutoDock 4 scoring function
            //// sUnbound_ls = call_ls( UnboundLocalSearchMethod, sUnbound_ext, ls_pop_size, &ligand );
            //// // sUnbound_ext = sUnbound_ls; // if you want to update sUnbound_ext to be sUnbound_ls...
            // --- Finished Local Search ---
            // end of Step 2 // }

            // Step 3 // {
            //
            // Restore the AutoDock 4 force field for docking
            // ----------------------------------------------
            //
            // Remember to turn on the use of the non-bond cutoff
            B_use_non_bond_cutoff = TRUE;
            //
            // Restore the setting for calculation of internal electrostatics to the saved value:
            B_calcIntElec = B_calcIntElec_saved;
            //
            // Restore the number of torsions in the state variables "sInit" and "sUnbound_ext"
            sInit.ntor = saved_sInit_ntor;
            sUnbound_ext.ntor = saved_sInit_ntor;
            //
            // Use the standard AutoDock energy tables to compute the internal energy
            // Use this value to set unbound_internal_FE
            evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom,
                            info, map, peratomE, nonbondlist, ad_energy_tables, Nnb,
			    Nnb_array, &group_energy,
                            B_calcIntElec, B_isGaussTorCon, B_isTorConstrained,
                            B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, sInit, ligand,
                            ignore_inter,
                            B_include_1_4_interactions, scale_1_4, scale_eintermol,
                            unbound_internal_FE, B_use_non_bond_cutoff, B_have_flexible_residues,
			    true_ligand_atoms, outlev, logFile);
                            //parameterArray, unbound_internal_FE, B_use_non_bond_cutoff, B_have_flexible_residues);
            // end of Step 3 // }

            // Step 4 // {
            //
            // Compute the energy of the unbound extended state
            // ------------------------------------------------
            //
            // Convert from unbound state to unbound coordinates
            cnv_state_to_coords( sUnbound_ext, vt, tlist, sUnbound_ext.ntor, crdpdb, crd, natom,
	     true_ligand_atoms, outlev, logFile);
            //
            // Calculate the unbound internal energy using the standard AutoDock energy function
            (void) eintcalPrint(nonbondlist, ad_energy_tables, crd, Nnb, Nnb_array, &group_energy,
	    B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge,
	    B_use_non_bond_cutoff, B_have_flexible_residues, natom, type, info->atom_type_name, outlev, logFile);
            //
            // eintcal() and eintcalPrint() set the values of group_energy[]
            unbound_ext_internal_FE = 
	      group_energy.intra_moving_moving_lig.total + 
	        group_energy.intra_moving_moving_rec.total;
            //
            pr(logFile, "\n\nThe internal energy of the unbound extended state was computed to be %+.3lf kcal/mol\n\n", unbound_ext_internal_FE);
            // end of Step 4 // }

            // Step 5 // {
            //
            // Decide whether to use extended or AutoDock state for unbound state
            // ------------------------------------------------------------------
            //
            if (unbound_ext_internal_FE > 0.0) {
                // Unbound extended state has an internal energy that is positive

                // Step 5.1 // {
                //
                // Repeat Step 1 with the standard AutoDock internal energy potentials
                //
                // Run a hybrid global-local search using the autodock energy tables
                // -----------------------------------------------------------------
                //
                //  *  Turn off the use of the non-bond cutoff
                //  *  Turn off internal electrostatics
                //  *  Turn off intermolecular energy calculations
                // TODO Need not translate or rotate the ligand in unbound searches
                //
                /*
                 *  Genetic Algorithm-Local search,  a.k.a.
                 *  Lamarckian Genetic Algorithm
                 */
                (void) fprintf( logFile, "\n\tBEGINNING COMPUTATION OF UNBOUND AUTODOCK STATE USING LGA\n");
                (void) fprintf( logFile,     "\t_________________________________________________________\n\n\n");
                //
                pr(logFile, "Date:\t");
                printdate( logFile, 2 );
                (void) fflush( logFile );
                //  Can get rid of the following line
                ((Genetic_Algorithm *)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor, outlev, logFile);
                //
                // Start Lamarckian GA run searching only torsions -- Unbound simulation
                // sUnbound_ad = call_glss_tors( GlobalSearchMethod, LocalSearchMethod,
                sUnbound_ad = call_glss( GlobalSearchMethod, LocalSearchMethod,
                                         sInit,
                                         num_evals_unbound, pop_size,
                                         outlev, logFile,
                                         output_pop_stats, &ligand, &evaluate,
                                         B_RandomTran0, B_RandomQuat0, B_RandomDihe0,
                                         info, FN_pop_file, end_of_branch );
                // State of best individual at end of GA-LS run, sUnbound_ad, is returned.
                // Finished Lamarckian GA run
                pr( logFile, "\nFinished Lamarckian Genetic Algorithm (LGA)\n");
                printdate( logFile, 1 );
                (void) fflush( logFile );
                pr(logFile, "\nTotal number of Energy Evaluations: %u\n", evaluate.evals() );
                pr(logFile, "Total number of Generations:        %u\n", ((Genetic_Algorithm *)GlobalSearchMethod)->num_generations());
                // Restore the number of torsions in the state variable "sUnbound_ad"
                sUnbound_ad.ntor = saved_sInit_ntor;
                // end of Step 5.1 // }

                // Step 5.2 // {
                //
                // Compute the energy of the unbound AutoDock state
                // ------------------------------------------------
                //
                // Convert from unbound state to unbound coordinates
                cnv_state_to_coords( sUnbound_ad, vt, tlist, sUnbound_ad.ntor, crdpdb, crd, natom,
		 true_ligand_atoms, outlev, logFile);
                //
                // Calculate the unbound internal energy using the standard AutoDock energy function
                (void) eintcalPrint(nonbondlist, ad_energy_tables, crd, Nnb, Nnb_array, &group_energy,
		B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, 
		B_use_non_bond_cutoff, B_have_flexible_residues, natom, type, info->atom_type_name, outlev, logFile);
                //
                // eintcal() and eintcalPrint() set the values of group_energy[]
                unbound_ad_internal_FE = 
	          group_energy.intra_moving_moving_lig.total + 
	            group_energy.intra_moving_moving_rec.total;
                //
                pr(logFile, "\n\nThe internal energy of the unbound AutoDock state was computed to be %+.3lf kcal/mol\n\n", unbound_ad_internal_FE);
                // end of Step 5.2 // }

                if (unbound_ad_internal_FE < unbound_ext_internal_FE) {
                    pr(logFile, "NOTE:   The AutoDock internal energy of the \"extended\" state was higher\nNOTE:   than that of the state obtained by searching using the AutoDock internal\nNOTE:   energy potentials.\nNOTE:   The unbound state was set to the AutoDock optimum state, not the \"extended\" state.\n\n");
                    unbound_internal_FE = unbound_ad_internal_FE;
                    sUnbound = sUnbound_ad;
                } else {
                    pr(logFile, "NOTE:   Although the AutoDock internal energy of the \"extended\" state was positive, it was lower\nNOTE:   than that of the state obtained by searching using the AutoDock internal\nNOTE:   energy potentials.\nNOTE:   The unbound state was set to the \"extended\" state.\n\n");
                    unbound_internal_FE = unbound_ext_internal_FE;
                    sUnbound = sUnbound_ext;
                }
            } else {
                // Unbound extended state has an internal energy that is negative
                unbound_internal_FE = unbound_ext_internal_FE;
                sUnbound = sUnbound_ext;
                pr(logFile, "NOTE:   The AutoDock internal energy of the \"extended\" state was negative.\n\nNOTE:   The unbound state was set to the \"extended\" state.\n\n");
            }
            //
            pr(logFile, "\n\nThe internal energy of the unbound state was set to %+.3lf kcal/mol\n\n", unbound_internal_FE);
            // end of Step 5 // }

            // Step 6 // {
            //
            // Convert from unbound state to unbound coordinates
            cnv_state_to_coords( sUnbound, vt, tlist, sUnbound.ntor, crdpdb, crd, natom,
	     true_ligand_atoms, outlev, logFile);
            // end of Step 6 // }

            // Step 7 // {
            //
            // Output the coordinates of the unbound state
            pr( logFile, "\n\n\tFINAL UNBOUND STATE\n" );
            pr( logFile,     "\t___________________\n\n\n" );
            //
            writePDBQT( -1, runseed[nconf],  FN_ligand, dock_param_fn, lig_center,
                        sUnbound, ntor, NULL, NULL, natom, atomstuff,
                        crd, peratomE,
                        charge, abs_charge, qsp_abs_charge,
                        ligand_is_inhibitor,
                        torsFreeEnergy,
                        vt, tlist, crdpdb, nonbondlist,
                        ad_energy_tables,
                        type, 
			Nnb, Nnb_array, &group_energy, true_ligand_atoms,
			B_calcIntElec,
                        map,
                        ignore_inter,
                        B_include_1_4_interactions, scale_1_4, parameterArray, unbound_internal_FE,
                        info, UNBOUND, PDBQT_record, B_use_non_bond_cutoff, B_have_flexible_residues, ad4_unbound_model,
			outlev, logFile);
            // end of Step 7 // }

            // Step 8 // {
            //
            // Remember to reset the energy evaluator back to computing the intermolecular energy between
            // the flexible and the rigid molecules.
            evaluate.compute_intermol_energy(TRUE);
            // end of Step 8 // }

        } else {
            pr(logFile, "NOTE:  AutoDock cannot compute the energy of the unbound state, since the ligand is rigid.\n\n");
            pr(logFile, "NOTE:  Use the \"unbound energy\" command to set the energy of the unbound state, if known from a previous calculation where the ligand was treated as flexible.\n\n");
            unbound_internal_FE = 0.0L;
            ad4_unbound_model = User;
            pr(logFile, "\n\nThe internal energy of the unbound state was set to %+.3lf kcal/mol\n\n", unbound_internal_FE);
        }

        pr( logFile, UnderLine );
        break;

/*____________________________________________________________________________*/

    case DPF_EPDB:
        /*
         *  epdb
         *
         *  Computes the energy of the ligand specified by the "move lig.pdbqt" command.
         *  Return the energy of the Small Molecule.
         *  FN_ligand must be in   PDBQT-format;
         *  flag can be:-
         *  0 = NEW, or   PDBQT-71, and
         *  1 = OLD, or   PDBQT-55 (old PDBq format).
         */
	{ // block for epdb locals:
	static EnergyComponent zeroEC; // const, always all zeros
	EnergyComponent totalE = zeroEC;
	Real emap_total = 0.; // does not include desolv
	Real desolv_total = 0.; 
	Real elec_total = 0.;
	Real charge_total = 0.;

	Real eintra=0;  // sum of intramolecular energy for the ligand plus that of the protein

        nfields = sscanf(line, "%*s %s", dummy_FN_ligand);
        if (nfields >= 1) {
            pr(logFile, "ERROR: \"epdb\" computes the energy of the ligand specified by the \"move lig.pdbqt\" command.\n");
            stop("it will not read in the PDBQT file specified on the \"epdb\" command line.");
        }

        exit_if_missing_elecmap_desolvmap_about("epdb");

        // warn if any atoms are outside the grid box
        for (i=0; i<natom; i++) {
	    Boole this_atom_outside;
	    this_atom_outside  = is_out_grid_info(crdorig[i][X], crdorig[i][Y], crdorig[i][Z]);
            if (this_atom_outside) {
                (void) sprintf( message, "%s: WARNING: Atom %d (%.3f, %.3f, %.3f) is outside the grid!\n", programname, i+1, crdorig[i][X], crdorig[i][Y], crdorig[i][Z] );
                print_2x( logFile, stderr, message );
            }
        }
        pr(logFile, "Number of \"true\" ligand atoms:  %d\n", true_ligand_atoms);
        //
        for (i=0;i<natom;i++) {
            if (ignore_inter[i] == 1) {
                pr(logFile, "Special Boundary Conditions:\n");
                pr(logFile, "____________________________\n\n");
                pr(logFile, "AutoDock will ignore the following atoms in the input PDBQT file \nin intermolecular energy calculations:\n");
                pr(logFile, "\n(This is because these residue atoms are at the boundary between \nflexible and rigid, and since they cannot move \nthey will not affect the total energy.)\n\n");
                break;
            }
        }
        for (i=0;i<natom;i++) {
            if (ignore_inter[i] == 1) {
                pr(logFile, "Atom number %d:  %s\n", i+1, atomstuff[i] );
            }
        }
        pr(logFile, "\n");

        sInit.ntor = ligand.S.ntor;

        // Calculate the internal energy
        if (ntor > 0) {
            eintra= eintcalPrint(nonbondlist, ad_energy_tables, crdorig, Nnb, Nnb_array, &group_energy,
	    B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, 
	    B_use_non_bond_cutoff, B_have_flexible_residues, natom, type, info->atom_type_name, outlev, logFile);
        }

        pr(logFile, "Unbound model to be used is %s.\n", report_parameter_library());
        // calculate the energy breakdown for the input coordinates, "crdorig"
        // Use 0.0 for the unbound internal free energy -- since this is a "single-point energy calculation"
        eb = calculateBindingEnergies( natom, ntor, 0.0 /*unbound_internal_FE*/, torsFreeEnergy, B_have_flexible_residues,
                                crdorig, charge, abs_charge, type, map, info,
                                ignore_inter, peratomE, &totalE,
                                nonbondlist, ad_energy_tables, Nnb, Nnb_array, &group_energy, true_ligand_atoms,
				B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, 
                                B_use_non_bond_cutoff, User /*ad4_unbound_model*/, outlev, logFile);

        pr(logFile, "\n\n\t\tIntermolecular Energy Analysis\n");
        pr(logFile,     "\t\t==============================\n\n");
        pr(logFile, "Atom Atom    Total      vdW+Hbond  Electrosta  Desolvation Partial          Coordinates         \n");
        pr(logFile, " Num Type    Energy      Energy    tic Energy     Energy    Charge      x         y         z    \n");
        pr(logFile, "____ ____  __________  __________  __________  __________ _______  ________  ________  ________\n");
        //          "1234  0123456789  0123456789  0123456789  1234567  12345678  12345678  12345678"

        for (int i = 0;  i < natom;  i++) {
            pr(logFile, "%4d  %-2s  %10.4f  %10.4f  %10.4f  %10.4f %8.4f  %8.4f  %8.4f  %8.4f\n",
	    i+1,  info->atom_type_name[type[i]],
	    peratomE[i].total, peratomE[i].vdW_Hb, peratomE[i].elec, peratomE[i].desolv, 
	    charge[i], crdorig[i][X], crdorig[i][Y], crdorig[i][Z]);

	    emap_total += peratomE[i].vdW_Hb;
	    elec_total += peratomE[i].elec;
	    desolv_total += peratomE[i].desolv;
            charge_total += charge[i];
        } /*i*/
        pr(logFile, "          __________  __________  __________  _______ ________\n");
        pr(logFile, "    Total %10.4lf  %10.4lf  %10.4lf %10.4f %8.4lf\n",        
	 (double)(emap_total + desolv_total + elec_total), 
	 (double)emap_total, 
	 (double)elec_total, 
	 (double)desolv_total, 
	 (double)charge_total);
        pr(logFile, "          __________  __________  __________  __________ _______\n");
        pr(logFile, "          vdW+Hbond+   vdW+Hbond  Electrosta  Desolvation Partial\n");
        pr(logFile, "          +Elec+Desolv   Energy    tic Energy    Energy    Charge\n");
        pr(logFile, "              Energy \n\n");

        pr(logFile, "Total Intermolecular Interaction Energy          = %+8.4lf kcal/mol\n", (double)eb.e_inter);
        pr(logFile, 
        "Total Intermolecular vdW + Hbond + desolv Energy = %+8.4lf kcal/mol\n",
	   (double) emap_total + desolv_total);
        pr(logFile, 
	 "Total Intermolecular Electrostatic Energy        = %+8.4lf kcal/mol\n", 
	   (double) elec_total );
        pr(logFile, "Total Intermolecular + Intramolecular Energy     = %+8.4lf kcal/mol\n", (double)eb.e_inter+eintra);
	pr(logFile, "\n\n");

	// see also writePDBQT.cc and analysis.cc for similar code:
        printEnergies( &eb, "epdb: USER    ", ligand_is_inhibitor, emap_total+desolv_total, elec_total, 
	 B_have_flexible_residues,  // next two terms are meaningful only if have flexible residues...
	 group_energy.inter_moving_moving.vdW_Hb + group_energy.inter_moving_moving.desolv,
	 group_energy.inter_moving_moving.elec,
	 User /*ad4_unbound_model*/,
	 outlev, logFile);
        pr(logFile, "\n");

	} // block for epdb locals:

        break;

/*____________________________________________________________________________*/

    case DPF_TERMINATION:
        /*
         *  ga_termination energy 0.1
         *  ga_termination evals 25000  // the best energy did not change in this time
         *  ga_termination time 120 s
         */
        /*
        (void) sscanf( line, "%*s %d", &i );
        */
        break;

/*____________________________________________________________________________*/

    case GA_CROSSOVER_MODE:
        /*
         * ga_crossover_mode OnePt
         * ga_crossover_mode TwoPt
         * ga_crossover_mode Uniform
         * ga_crossover_mode Arithmetic
         *
         * Xover_Mode c_mode = OnePt;  //  can be: OnePt, TwoPt, Uniform or Arithmetic
         */
	c_mode_str[0]='\0';
        get1arg( line, "%*s %s", c_mode_str, "GA_CROSSOVER_MODE" );
        if (streq(c_mode_str, "onept")) {
            c_mode = OnePt;
            pr(logFile, "One-point crossover will be used in GA and LGA searches.\n");
        } else if (streq(c_mode_str, "twopt")) {
            c_mode = TwoPt;
            pr(logFile, "Two-point crossover will be used in GA and LGA searches.\n");
        } else if (streq(c_mode_str, "uniform")) {
            c_mode = Uniform;
            pr(logFile, "Uniform crossover will be used in GA and LGA searches.\n");
        } else if (streq(c_mode_str, "arithmetic")) {
            c_mode = Arithmetic;
            pr(logFile, "Arithmetic crossover will be used in GA and LGA searches.\n");
        } else if (streq(c_mode_str, "branch")) {
            c_mode = Branch;
            pr(logFile, "Branch crossover will be used in GA and LGA searches.\n");
        } else stop("unrecognized mode in GA_CROSSOVER_MODE line");
        break;

/*____________________________________________________________________________*/

    case GA_TOURNAMENT_SELECTION:
        /*
         * ga_tournament_selection
         */
        
	// M Pique - does not appear to work so disabling for now (Oct 2009)
        //s_mode = Tournament;
        //pr(logFile, "Tournament selection will be used in GA and LGA searches.\n");
         prStr( error_message, "%s:  ERROR! Tournament selection is not yet implemented!\n", programname);
         pr_2x( logFile, stderr, error_message );
	 stop(error_message);
         break;

/*____________________________________________________________________________*/

    case GA_LINEAR_RANKING_SELECTION:
        /*
         * ga_linear_ranking_selection [probability_ratio Real]
         */
        
	// look for optional ratio value
        nfields = sscanf( line, "%*s %s", c_mode_str );
        if (nfields>0&&streq(c_mode_str, "probability_ratio")){ 
            nfields = sscanf( line, "%*s %*s " FDFMT, &linear_ranking_selection_probability_ratio);
	    if(nfields!=1) stop("syntax error in GA_LINEAR_RANKING_SELECTION line");
        }

        s_mode = LinearRanking;
        pr(logFile, "Linear ranking selection will be used in GA and LGA searches.\n");
        break;

/*____________________________________________________________________________*/

    case GA_PROPORTIONAL_SELECTION:
        /*
         * ga_proportional_selection 
         */
        
        s_mode = Proportional;
        pr(logFile, "Proportional selection will be used in GA and LGA searches.\n");
        break;

/*____________________________________________________________________________*/

    case GA_BOLTZMAN_SELECTION:
        /*
         * ga_boltzman_selection 
         */
        
        prStr( error_message, "%s:  ERROR! Boltzman selection is not yet implemented! \n", programname);
        pr_2x( logFile, stderr, error_message );
	stop(error_message);
        break;

/*____________________________________________________________________________*/


    case DPF_POPFILE:
        /*
         *  output_pop_file
         *
         *  Used to write out the population to a file at the end of
         *  every GA.
         */
        get1arg( line, "%*s %s", FN_pop_file, "OUTPUT_POP_FILE");
        pr( logFile, "The population will be written to the file \"%s\" at the end of every generation.\n", FN_pop_file);
        break;

 /*____________________________________________________________________________*/

     case DPF_CONFSAMPLER:
        /*
         * confsampler
         *
         * Scan a region around conformations saved in sHist array.
         *
         */

        nfields = sscanf( line, "%*s %s %d", confsampler_type, &confsampler_samples);
	if(nfields<1) stop("syntax error in CONFSAMPLER line");
        pr( logFile, "Scanning local regions around each docked conformation.\n");

        exit_if_missing_elecmap_desolvmap_about("confsampler");

        if (streq(confsampler_type, "systematic")) {
            systematic_conformation_sampler(sHist, nconf, vt, crdpdb, tlist,
	     lig_center, natom, type, info, true_ligand_atoms, &evaluate, outlev, logFile);
        } 
	else if (streq(confsampler_type, "random")) {
	    if(nfields<2) stop("syntax error in CONFSAMPLER RANDOM line");
            random_conformation_sampler(sHist, nconf, confsampler_samples, vt, crdpdb, tlist,
	    lig_center, natom, type, info, true_ligand_atoms, &evaluate, outlev, logFile);
        }
        else stop("unrecognized mode in in CONFSAMPLER line");
        break;

/*____________________________________________________________________________*/

    case DPF_COPYRIGHT:
        /*
         * 'copyright' to show the AutoDock copyright notice
         */
        show_copyright(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_WARRANTY:
        /*
         * 'warranty' to show the AutoDock warranty
         */
        show_warranty(logFile);
        break;

//______________________________________________________________________________
/* illegal or unhandled keyword - fatal error as of May 2010 */

    case DPF_UNKNOWN:
    default:
        /*
         *
         */
         prStr(error_message, "unrecognized DPF keyword in \"%s\"", line);
         stop(error_message);
	break;

//______________________________________________________________________________

    } /* switch( dpf_keyword ) */

} /* while PARSING-DPF parFile */

/* __________________________________________________________________________
**
** Close the docking parameter file...
** __________________________________________________________________________
*/
pr( logFile, ">>> Closing the docking parameter file (DPF)...\n" );
(void) fclose( parFile );


//______________________________________________________________________________
/*
** Print the time and date when the docking has finished...
*/

pr( logFile, "This docking finished at:\t\t\t" );
printdate( logFile, 1 );
pr( logFile, "\n" );

success( hostnm, jobStart, tms_jobStart, logFile);

 if(write_stateFile){
   fprintf(stateFile,"</autodock>\n");
   (void) fclose( stateFile );
 }
(void) fclose( logFile );

// delete arrays
delete []nonbondlist;

//________________________________________________________________________________
/*
** End of Boinc
*/
#ifdef BOINCCOMPOUND
 boinc_fraction_done(1.);
#endif
#ifdef BOINC
    boinc_finish(0);       /* should not return */
#endif

return EXIT_SUCCESS;

} /* END OF PROGRAM */

/* AutoDock main private utility functions
*/
static void exit_if_missing_elecmap_desolvmap_about(string keyword)
{

    char error_message[LINE_LEN+100];

    if (!B_found_elecmap) {
         prStr(error_message, "%s:  %s command: no \"elecmap\" command has been specified!\n", programname, keyword.c_str());
         stop(error_message);
         exit(EXIT_FAILURE);
    } else if (!B_found_desolvmap) {
         prStr(error_message, "%s:  %s command: no \"desolvmap\" command has been specified!\n", programname, keyword.c_str());
         stop(error_message);
         exit(EXIT_FAILURE);
    }
}

static int getoutlev(char *line, int *outlev) {
 /* set *outlev either numerically or symbolically. return 1 if OK, 0 for fail
  * Allow either integers or symbolic values to appear first or second in line
  *  e.g.   "outlev 1", "1", "outlev adt", "ADT"    with case not significant
  */
	char s[LINE_LEN];

	if(1==sscanf(line, "%*s %d", outlev)||1==sscanf(line, "%d", outlev)) return 1;
	// if not integer, look for symbolic outlev, see constants.h
	if(1!=sscanf(line, "%*s %s", s) && 1!= sscanf(line, "%s", s)) return 0;
	for(int i=0;strlen(outlev_lookup[i].key)>0;i++) { 
	  if(streq(s,outlev_lookup[i].key)) {
	    *outlev=outlev_lookup[i].value;
	    return 1;
	  }
	}
	return 0;
}
static void set_seeds( FourByteLong seed[2], char seedIsSet[2], FourByteLong runseed[MAX_RUNS][2], const int outlev, FILE *logFile ) {
	    time_t time_seed;
	 if ( (!seedIsSet[0]) && (sizeof seed[0])!=4) 
	   pr(logFile, "Sizeof FourByteLong = %d, not 4.\n",
	   (int) (sizeof seed[0])); // DEBUG 2012-07 M Pique
	 if( ! seedIsSet[0] ) seed[0] = (FourByteLong) processid();
	 if( ! seedIsSet[1] ) seed[1] = (FourByteLong) time( &time_seed );
	 seedIsSet[0]  = seedIsSet[1] = 'D';

	 /* set seeds now for all possible runs - so each run's results will be
          * independent of other runs, for parallelization 
	  * Note that the first run will use the DPF-specified seed, if any,
	  * possibly as modified by population setup work, see comment above.
	  * The setall() assures that every random number generator (thread)
	  * will have these seeds.
	  */
         setall(seed[0], seed[1]);  // see com.cc 
	 runseed[0][0] = seed[0];
	 runseed[0][1] = seed[1];
         for(int j=1;j<MAX_RUNS;j++) {
		runseed[j][0] = ignlgi();
		runseed[j][1] = ignlgi();
		}
		
	 /* make sure RNG #0 is set back to DPF-specified seed for compatibility with existing tests */
	 (void) gscgn(0);
	 setsd_t( seed[0], seed[1],0 );

         if(outlev>=LOGMIN) pr(logFile,
	  "Random number generator was seeded with values " FBL_FMT ", " FBL_FMT ".\n",
	  seed[0], seed[1]);
	 }

static int processid() {
#ifdef HAVE_GETPID
                    return getpid();
#elif HAVE_GETPROCESSID
                    return GetProcessId(); // Windows WIN32
#else
		    pr(logFile, "Cannot determine process id for random number seed, using dummy '12345'");
		    return 12345;
#endif
}

#ifdef BOINC
/*  Dummy graphics API entry points.
 *  This app does not do graphics, but it still must provide these callbacks.
 */

void app_graphics_render(int xs, int ys, double time_of_day) {}
void app_graphics_reread_prefs(){}
void boinc_app_mouse_move(int x, int y, bool left, bool middle, bool right ){}
void boinc_app_mouse_button(int x, int y, int which, bool is_down){}
void boinc_app_key_press(int wParam, int lParam){}
void boinc_app_key_release(int wParam, int lParam){}
#endif

/* EOF */
