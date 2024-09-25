#! /usr/bin/env python
#
# $Id: test_autodock4.py,v 1.74 2014/07/16 23:27:27 mp Exp $
#

"""
Unit Tests for AutoDock 4.
"""

#______________________________________________________________________________

import sys
import os
import unittest
import getopt
import subprocess
from DlgParser import DlgParser

#______________________________________________________________________________
#
# Global variables

autodock_executable = "../autodock4" # where the AutoDock executable resides
dpf_directory = '.' # where the input DPF files reside
test_output_directory = '.' # where the DLG files will be written
fnull = open(os.devnull, "w")

try:
    opts, argv = getopt.getopt(sys.argv[1:], "d:e:o:",
    ["dpf-directory=","executable=","test-output-directory="]) 
except getopt.GetoptError, v:
    usage() 
    sys.exit(2)

for o,a in opts:
    if o in ("-d", "--dpf-directory"):
        dpf_directory = a
    if o in ("-e", "--executable"):
        autodock_executable = a
    if o in ("-o","--test-output-directory"):
        test_output_directory = a


#______________________________________________________________________________

def usage():
    """Print out the usage of this command."""
    print """Usage:  python test_autodock4.py [-d <string>] [-e <string>] [-o <string>]

where:
    -d, --dpf-directory
        specifies the directory containing the DPFs to be tested;
        this flag is optional; default is '.'
    -e, --executable
        specifies the path to the AutoDock executable to be tested;
        this flag is optional; default is '../autodock4'
    -o, --test-output-directory
        specifies the directory where the output DLGs will be written;
        this flag is optional; default is '.'

NOTE:  these may be relative to the directory where this script was invoked.
"""

#______________________________________________________________________________

def run_AutoDock( dpf_filename, dlg_filename ):
    """Launch AutoDock, using the specified AutoDock executable and DPF,
    create the specified DLG. Silently discard all output from standard
    output and standard error, so users are not bothered by the expected
    and routine AutoDock error messages from many of the tests."""
    dpf = dpf_directory + os.sep + dpf_filename
    dlg = test_output_directory + os.sep + dlg_filename
    rm( dlg )
    command =   [autodock_executable, '-p', dpf, '-l', dlg ] 
    print '\nRunning ' + autodock_executable + ' using DPF "'+dpf+'", saving results in "'+dlg+'":'
    try:
        rc = subprocess.call( command, stdout = fnull, stderr = fnull )
	#print 'autodock returned ', rc  # DEBUG
        return find_success_in_DLG( dlg_filename )
    except OSError,e:
        print "\nUnable to run " + autodock_executable + " :", e
        return False

#______________________________________________________________________________

def parse_energy_from_DLG( dlg_filename, energy_list):
    """Parse the AutoDock DLG, and return the intermolecular and internal
    energies as a tuple."""
    parser = DlgParser()
    dlg = test_output_directory + os.sep + dlg_filename
    parser.parse( dlg )
    docked = parser.clist[0]  #dictionary of results
    result = []
    for energy_type in energy_list:
        newVal = docked.get(energy_type, 'ERROR')
        print energy_type, ' is now ', newVal
        result.append(docked.get(energy_type, 'ERROR'))
    #intermol_energy = docked['intermol_energy']  #-6.17
    #internal_energy = docked['total_internal']  # -1.58
    #print "docked[binding_energy]=", docked['binding_energy']
    #print "docked[electrostatic_energy]=", docked['electrostatic_energy']
    #print "docked[intermol_energy]=", docked['intermol_energy']
    #print "docked[total_internal]=", docked['total_internal']
    #unbound_energy = docked['unbound_energy']
    #print "unbound_energy=", unbound_energy
    #return ( intermol_energy, internal_energy )
    return result

#______________________________________________________________________________

def find_success_in_DLG( dlg_filename ):
    """Open the AutoDock DLG, and look for the string "Successful Completion"
    in the last 10 lines of the file."""
    dlg = test_output_directory + os.sep + dlg_filename
    try:
        fptr = open( dlg )
        lines = fptr.readlines()
        fptr.close()
        success = False
        for l in lines[-10:]:
            if l.find( "Successful Completion" ) > -1:
                success = True
        return success
    except:
        return False

#______________________________________________________________________________

def rm( filename ):
    """Remove single file """
    # no regular expression wildcards, e.g. * or ?, allowed
    if os.access(filename, os.F_OK) :
       try:
          os.remove( filename )
       except OSError,e:
         print 'unable to remove '+filename


#______________________________________________________________________________

class AutoDock_base_test( unittest.TestCase ):
    """Base Class for AutoDock testing."""
    """ do not instantiate this for a test, use simple_test or other subclass instead """
    dpf_stem = "BaseClass"
    computed = False
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        self.dlg_filename = "test_" + self.dpf_stem + ".dlg"
        self.computed = run_AutoDock( self.dpf_stem + ".dpf", self.dlg_filename )

    #def test_dlg_exists( self ):
    #    """Check that run finished and a new DLG has been computed."""
    #    # Check that run finished and a new DLG has been computed.
    #    if (self.expected_outcome == True ):
    #        print "Testing that DLG exists and AutoDock successfully completed."
    #    else:
    #        print "Testing that DLG exists and AutoDock did not complete."
    #    self.assertEqual( self.computed, self.expected_outcome )


#______________________________________________________________________________

class AutoDock_simple_test( unittest.TestCase ):
    """Base Class for AutoDock testing."""
    dpf_stem = "BaseClass"
    computed = False
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        self.dlg_filename = "test_" + self.dpf_stem + ".dlg"
        self.computed = run_AutoDock( self.dpf_stem + ".dpf", self.dlg_filename )

    def test_dlg_exists( self ):
        """Check that run finished and a new DLG has been computed."""
        # Check that run finished and a new DLG has been computed.
        if (self.expected_outcome == True ):
            print "Testing that DLG exists and AutoDock successfully completed."
        else:
            print "Testing that DLG exists and AutoDock did not complete."
        self.assertEqual( self.computed, self.expected_outcome )
#______________________________________________________________________________

class AutoDock4_1pgp_no_extension( AutoDock_simple_test ):
    """Test that autodock4 stops early if .dpf extension is missing
    keywords are specified."""
    dpf_stem = "1pgp_no_extension"
    print "in 1pgp_no_extension"
    expected_outcome = False # True means Successful Completion!
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        self.dlg_filename = "test_" + self.dpf_stem + ".dlg"
        self.computed = run_AutoDock( self.dpf_stem , self.dlg_filename )

#______________________________________________________________________________

class AutoDock4_1pgp_wrong_extension( AutoDock_simple_test ):
    """Test that autodock4 stops early if extension is not '.dpf'
    keywords are specified."""
    dpf_stem = "1pgp.fpd"
    print "in 1pgp_wrong_extension"
    expected_outcome = False # True means Successful Completion!
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        self.dlg_filename = "test_" + self.dpf_stem + ".dlg"
        self.computed = run_AutoDock( self.dpf_stem , self.dlg_filename )

#______________________________________________________________________________


class AutoDock4_1pgp_two_extensions( AutoDock_simple_test ):
    """Test that autodock4 stops early if dpf name includes two .dpf
    keywords are specified."""
    dpf_stem = "1pgp.dpf"
    print "in 1pgp_two_extensions"
    expected_outcome = False # True means Successful Completion!
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        self.dlg_filename = "test_" + self.dpf_stem + ".dlg"
        self.computed = run_AutoDock( self.dpf_stem +'.dpf', self.dlg_filename )


#______________________________________________________________________________

class AutoDock4_1pgp_ligand_types_map_mismatch( AutoDock_simple_test ):
    """Test that autodock4 stops early if number of maps do not equal number
    of ligand types"""
    dpf_stem = "1pgp_ligand_types_map_mismatch"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________
class AutoDock4_1pgp_map_truncated( AutoDock_simple_test ):
    """Test that autodock4 stops early if map data is incomplete"""
    dpf_stem = "1pgp_map_truncated"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_illegal_keyword_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if it finds an illegal keyword 
    in dpf """
    dpf_stem = "1pgp_illegal_keyword"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_no_elecmap_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if no "elecmap" keyword is specified."""
    dpf_stem = "1pgp_no_elecmap"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_no_desolvmap_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if no "desolvmap" keyword is specified."""
    dpf_stem = "1pgp_no_desolvmap"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_no_elec_desolv_maps_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if no elecmap and no desolvmap 
    keywords are specified."""
    dpf_stem = "1pgp_no_elec_desolv_maps"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_too_many_torsions( AutoDock_simple_test ):
    """Test that autodock4 stops early if too many torsions 
    are specified. (current limit is 32)"""
    dpf_stem = "1pgp_too_many_torsions"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_just_right_number_torsions( AutoDock_simple_test ):
    """Test that autodock4 completes with current limit of number of  torsions 
    are specified. (current limit is 32)"""
    dpf_stem = "1pgp_just_right_number_torsions"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________


class AutoDock4_1pgp_too_many_ligand_types_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if too many ligand types 
    are specified. (current limit is 14)"""
    dpf_stem = "1pgp_too_many_ligand_types"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_two_ligands_test( AutoDock_simple_test ):
    """Test that autodock4 can run dpf specifying two ligands."""
    dpf_stem = "1pgp_two_ligands"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_two_mapsets_test( AutoDock_simple_test ):
    """Test that autodock4 can run dpf specifying two sets of maps and one ligand."""
    dpf_stem = "1pgp_two_mapsets"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_set_illegal_test( AutoDock_simple_test ):
    """Test that autodock 4.1 works when unbound is set to 'foo' in the DPF."""
    dpf_stem = "1pgp_unbound_set_illegal"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_model_illegal_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if it finds an illegal unbound_model 
    in dpf """
    dpf_stem = "1pgp_unbound_model_illegal"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_mixed_adt_analysis_test( AutoDock_simple_test ):
    """Test that autodock 4 with outlevel ADT runs mixed between search methods followed by cluster analysis"""
    dpf_stem = "1pgp_mixed_adt_analysis"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_mixed_basic_analysis_test( AutoDock_simple_test ):
    """Test that autodock 4 with outlevel BASIC runs mixed between search methods followed by cluster analysis"""
    dpf_stem = "1pgp_mixed_basic_analysis"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_ga_run_maxruns_test( AutoDock_simple_test ):
    """Test that autodock 4 runs with max number of runs"""
    dpf_stem = "1pgp_ga_run_maxruns"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_ga_run_maxruns_analysis_test( AutoDock_simple_test ):
    """Test that autodock 4 runs with max number of runs followed by cluster analysis"""
    dpf_stem = "1pgp_ga_run_maxruns_analysis"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________
class AutoDock4_1pgp_mixed_maxruns_test( AutoDock_simple_test ):
    """Test that autodock 4 runs with max number of runs mixed between search methods"""
    dpf_stem = "1pgp_mixed_maxruns"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_mixed_maxruns_analysis_test( AutoDock_simple_test ):
    """Test that autodock 4 runs with max number of runs  mixed between search methods followed by cluster analysis"""
    dpf_stem = "1pgp_mixed_maxruns_analysis"
    expected_outcome = True # True means Successful Completion!


#______________________________________________________________________________

class AutoDock4_1pgp_mixed_maxruns_adt_analysis_test( AutoDock_simple_test ):
    """Test that autodock 4 runs with max number of runs  mixed between search methods followed by cluster analysis with outlev 1 ADT"""
    dpf_stem = "1pgp_mixed_maxruns_adt_analysis"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_mixed_maxruns_basic_analysis_test( AutoDock_simple_test ):
    """Test that autodock 4 runs with max number of runs  mixed between search methods followed by cluster analysis with outlev 0 BASIC"""
    dpf_stem = "1pgp_mixed_maxruns_basic_analysis"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_overmaxruns_test( AutoDock_simple_test ):
    """Test that autodock 4 stops with over max number of runs  mixed between search methods"""
    dpf_stem = "1pgp_overmaxruns"
    expected_outcome = False # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_intnbpreps_test( AutoDock_simple_test ):
    """Test that autodock 4 handles intnbp_r_eps OK """
    dpf_stem = "1pgp_intnbpreps"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_intnbpreps_toofew_test( AutoDock_simple_test ):
    """Test that autodock 4 handles intnbp_r_eps with too few tokens """
    dpf_stem = "1pgp_intnbpreps_toofew"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_intnbpreps_illegal_atom_test( AutoDock_simple_test ):
    """Test that autodock 4 handles intnbp_r_eps with bad atom name """
    dpf_stem = "1pgp_intnbpreps_illegal_atom"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_intnbcoeffs_test( AutoDock_simple_test ):
    """Test that autodock 4 handles intnbp_coeffs OK """
    dpf_stem = "1pgp_intnbcoeffs"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_ga_select_tournament_test( AutoDock_simple_test ):
    """Test that autodock 4.2 stops when ga_select_tournament is set in the DPF."""
    dpf_stem = "1pgp_ga_select_tournament"
    expected_outcome = False # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_ga_select_linear_ranking_test( AutoDock_simple_test ):
    """Test that autodock 4.2 works when ga_select_linear_ranking is set in the DPF."""
    dpf_stem = "1pgp_ga_select_linear_ranking"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

#print "docked[binding_energy]=", docked['binding_energy']
class AutoDock4_1pgp_ligrand_ga_only_test( AutoDock_simple_test ):
    """Test that autodock 4.2 works when ligand is randomized within autodock as set in the DPF."""
    dpf_stem = "1pgp_ligrand_ga_only"
    #print "in new ligrand_ga_only test"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_ga_only_test( AutoDock_simple_test ):
    """Test that autodock 4.2 works when ga_only is set in the DPF."""
    dpf_stem = "1pgp_ga_only"
    #print "in new ga_only test"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_about_only_test( AutoDock_simple_test ):
    """Test that autodock 4.2 works when about is set in the DPF but tran0,
    dihe0 and quat0 are missing."""
    dpf_stem = "1pgp_about_only"
    #print "in new ga_only test"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_multipleabouts_test( AutoDock_simple_test ):
    """Test that autodock 4.2 works when about appears multiple times"""
    dpf_stem = "1pgp_multipleabouts"
    #print "in new ga_only test"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________


class AutoDock4_1pgp_simanneal_small( AutoDock_simple_test ):
    """Test that autodock 4.2 simanneal works (at all). """
    dpf_stem = "1pgp_simanneal_small"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_simanneal_small_geometric_schedule( AutoDock_simple_test ):
    """Test that autodock 4.2 simanneal alternate cooling schedule """
    dpf_stem = "1pgp_simanneal_small_geometric_schedule"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_about_auto_simanneal( AutoDock_simple_test ):
    """Test that autodock 4.2 simanneal works when about, tran0, dihe0, and quat0 are missing."""
    dpf_stem = "1pgp_about_auto_simanneal"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_about_auto_ga_only( AutoDock_simple_test ):
    """Test that autodock 4.2 ga_only works when about, tran0, dihe0, and quat0 are missing."""
    dpf_stem = "1pgp_about_auto_ga_only"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_about_auto_gals( AutoDock_simple_test ):
    """Test that autodock 4.2 ga_ls works when about, tran0, dihe0, and quat0 are missing."""
    dpf_stem = "1pgp_about_auto_gals"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_about_auto_local_only( AutoDock_simple_test ):
    """Test that autodock 4.2 local_only works when about, tran0, dihe0, and quat0 are missing."""
    dpf_stem = "1pgp_about_auto_local_only"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_flexres_test( AutoDock_simple_test ):
    """Test that autodock 4.2 works when flexres is set in the DPF."""
    dpf_stem = "1pgp_flexres"
    #print "in new flexres test"
        #expected_binding_energy =  -4.22
        #expected_intermol_energy = -5.93
        #expected_internal_energy = -2.18
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_rmsmode_heavy_atoms_only_noH_test( AutoDock_simple_test ):
    """Test that autodock4 works if noHs when using new keyword rmsmode heavy_atoms_only."""
    dpf_stem = "1pgp_rmsmode_heavy_atoms_only_noH"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_rmsmode_all_heavy_atom_pairs_only_test( AutoDock_simple_test ):
    """Test that autodock4 works if rms calc uses all pairs of heavy_atoms.
       Set by including these TWO rmsmode cmds in dpf:  
        rmsmode unique_pair       #sets B_unique_pair_flag to TRUE (default is FALSE)
        rmsmode heavy_atoms_only  #sets B_rms_heavy_atoms_only to TRUE (default is FALSE)
        """
    dpf_stem = "1pgp_rmsmode_all_heavy_atom_pairs_only"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_rmsmode_unique_heavy_atom_pairs_only_test( AutoDock_simple_test ):
    """Test that autodock4 works if rms calc uses only unique pairs of heavy_atoms specified 
          with two rmsmode cmds: unique_pair heavy_atoms_only."""
    dpf_stem = "1pgp_rmsmode_unique_heavy_atom_pairs_only"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock_test( AutoDock_base_test ):
    """Class for AutoDock testing."""

    def test_dlg_exists_and_test_energy( self ):
        """Check that run finished and a new DLG has been computed.
        Also check the final energy is the expected value."""
        # Check that run finished and a new DLG has been computed.
        if (self.expected_outcome == True ):
            print "Testing that DLG exists and AutoDock successfully completed."
        else:
            print "Testing that DLG exists and AutoDock did not complete."
        self.assertEqual( self.computed, self.expected_outcome )
        # Check the final energy is expected value.
        # These values are for the quick GALS search in 1pgp.dpf and relatives
        expected_intermol_energy = -6.17  # -6.44 for Real==float, -6.17 for Real==double
        expected_internal_energy = -3.23  # -3.28 for Real==float, -3.23 for Real==double
        (intermol_energy, internal_energy) = parse_energy_from_DLG( self.dlg_filename, ['intermol_energy','total_internal'] )
        print "Testing that intermolecular energy = %.2f kcal/mol." % (expected_intermol_energy,)
        self.assertEqual( round(intermol_energy,6), round(expected_intermol_energy,6))
        print "Testing that internal energy = %.2f kcal/mol." % (expected_internal_energy,)
        self.assertEqual( round(internal_energy,6), round(expected_internal_energy,6))
#______________________________________________________________________________

class AutoDock4_1pgp_test( AutoDock_test ):
    """Test that autodock4 executes using an extremely short run."""
    dpf_stem = "1pgp"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

# next test not installed yet as I do not have a way I know of to pass
# different expected values to the superclass - M Pique  July 2012
class AutoDock4_1pgp_intelec_off_test( AutoDock_test ):
    """Test autodock4.3 internal energy off option yields 4.2 value."""
    dpf_stem = "1pgp_intelec_off"
    expected_intermol = -7.72
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_timepid_test( AutoDock_simple_test ):
    """Test that autodock4 executes using seed set by time and process-id"""
    dpf_stem = "1pgp_timepid"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_seed01_test( AutoDock_simple_test ):
    """Test that autodock4 stops using seeds  set to 0 1"""
    dpf_stem = "1pgp_seed01"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_seed37_test( AutoDock_simple_test ):
    """Test that autodock4 executes using seeds  set to 3 7"""
    dpf_stem = "1pgp_seed37"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_seedinttime_test( AutoDock_simple_test ):
    """Test that autodock4 executes using seeds  set to integer and time"""
    dpf_stem = "1pgp_seedinttime"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_epdb_flexres_intelec_test( AutoDock_simple_test ):
    """Test autodock 4.3 epdb output flexres with intelec on."""
    dpf_stem = "1pgp_epdb_flexres_intelec"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_epdb_flexres_nointelec_test( AutoDock_simple_test ):
    """Test autodock 4.3 epdb output flexres with intelec off."""
    dpf_stem = "1pgp_epdb_flexres_nointelec"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_epdb_etables_test( AutoDock_simple_test ):
    """Test autodock 4.3 high outlev detail for epdb and energy tables."""
    dpf_stem = "1pgp_epdb_etables"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_numoutlev_test( AutoDock_simple_test ):
    """Test autodock 4.3 numeric outlev detail settings."""
    dpf_stem = "1pgp_numoutlev"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_no_outlev_test( AutoDock_simple_test ):
    """Test autodock 4.3 default outlev detail setting."""
    dpf_stem = "1pgp_no_outlev"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_illegal_outlev_test( AutoDock_simple_test ):
    """Test autodock 4.3 outlev with no argument setting."""
    dpf_stem = "1pgp_illegal_outlev"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_symoutlev_test( AutoDock_simple_test ):
    """Test autodock 4.3 symbolic and numeric outlev detail settings."""
    dpf_stem = "1pgp_symoutlev"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_illegal_symoutlev_test( AutoDock_simple_test ):
    """Test autodock 4.3 symbolic outlev incorrect setting."""
    dpf_stem = "1pgp_illegal_symoutlev"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_gals_set_sw1_test( AutoDock_simple_test ):
    """Test that autodock4 gals_run succeeds with older sw1 rather than psw1."""
    dpf_stem = "1pgp_gals_set_sw1"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_smaller_test( AutoDock_simple_test ):
    expected_intermol_energy = -6.17  
    expected_internal_energy = -3.23
    """Test that autodock4 executes using fewer parameters and an extremely short run."""
    dpf_stem = "1pgp_smaller"
    #expected_intermol_energy = -6.17  
    #expected_internal_energy = -3.23
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_gals_use_defaults_test( AutoDock_test ):
    """Test that autodock4 executes using default gals parameters."""
    dpf_stem = "1pgp_gals_use_defaults"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________
class AutoDock4_1pgp_gals_use_defaults_test( AutoDock_test ):
    """Test that autodock4 executes using default gals parameters."""
    dpf_stem = "1pgp_gals_use_defaults"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_gals_10runs_test( AutoDock_simple_test ):
    """Test that autodock4 executes 'real world' gals parameters."""
    dpf_stem = "1pgp_gals_10runs"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_simanneal_10runs_test( AutoDock_simple_test ):
    """Test that autodock4 executes using 'real world' SA parameters."""
    dpf_stem = "1pgp_simanneal_10runs"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_lsonly_10runs_test( AutoDock_simple_test ):
    """Test that autodock4 executes using 'real world' local search parameters."""
    dpf_stem = "1pgp_lsonly_10runs"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________
class AutoDock4_1pgp_no_set_ga_test( AutoDock_simple_test ):
    """Test that autodock4 ga_run fails when set_ga command is omitted."""
    dpf_stem = "1pgp_no_set_ga"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_no_set_psw1_test( AutoDock_simple_test ):
    """Test that autodock4 gals_run fails when set_psw1 command is omitted."""
    dpf_stem = "1pgp_no_set_psw1"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_no_parameter_file_test( AutoDock_test ):
    """Test that autodock4 works using default parameter library."""
    dpf_stem = "1pgp_no_parameter_file"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_output_pop_file_test( AutoDock_simple_test ):
    """Test that autodock4 executes while writing out popfile"""
    dpf_stem = "1pgp_outpopfil"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_rmsmode_heavy_atoms_only_test( AutoDock_test ):
    """Test that autodock4 works using new keyword rmsmode heavy_atoms_only."""
    dpf_stem = "1pgp_rmsmode_heavy_atoms_only"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________
class AutoDock4_energy_test( AutoDock_base_test ):
    """Class for AutoDock testing free energy."""
    expected_binding_energy = None

    def test_dlg_exists_and_test_energy( self):
        """Check that run finished and a new DLG has been computed.
        Also check the final energy is the expected value."""
        # Check that run finished and a new DLG has been computed.
        if (self.expected_outcome == True ):
            print "Testing that DLG exists and AutoDock successfully completed."
        else:
            print "Testing that DLG exists and AutoDock did not complete."
        self.assertEqual( self.computed, self.expected_outcome )
        # Check the final energy is expected value.
        #expected_unbound_energy = -3.23
        #expected_binding_energy = +843.59     -5.87
        (binding_energy) = parse_energy_from_DLG( self.dlg_filename, ['binding_energy'])[0]
        print "Testing that binding energy = %.2f kcal/mol." % (self.expected_binding_energy,)
        print "binding_energy=", binding_energy
        self.assertEqual( round(binding_energy,6), round(self.expected_binding_energy,6))
#______________________________________________________________________________


class AutoDock4_unbound_test( AutoDock_base_test ):
    """Class for AutoDock testing unbound energy."""
    expected_unbound_energy = None

    def test_dlg_exists_and_test_energy( self):
        """Check that run finished and a new DLG has been computed.
        Also check the final energy is the expected value."""
        # Check that run finished and a new DLG has been computed.
        if (self.expected_outcome == True ):
            print "Testing that DLG exists and AutoDock successfully completed."
        else:
            print "Testing that DLG exists and AutoDock did not complete."
        self.assertEqual( self.computed, self.expected_outcome )
        # Check the final energy is expected value.
        #expected_unbound_energy = -3.23
        (unbound_energy) = parse_energy_from_DLG( self.dlg_filename, ['unbound_energy'])[0]
        print "Testing that unbound energy = %.2f kcal/mol." % (self.expected_unbound_energy,)
        print "unbound_energy=", unbound_energy
        self.assertEqual( round(unbound_energy,6), round(self.expected_unbound_energy,6))
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_default_test( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound is NOT set in the DPF."""
    dpf_stem = "1pgp_unbound_default"
    expected_unbound_energy = -3.23
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_model_extended( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound_model is set to extended."""
    dpf_stem = "1pgp_unbound_model_extended"
    expected_unbound_energy = -0.75
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_compute_unbound_extended( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound_model is set to extended."""
    dpf_stem = "1pgp_unbound_compute_unbound_extended"
    expected_unbound_energy = -0.75
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_model_value( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound_model is set to a value in the DPF."""
    dpf_stem = "1pgp_unbound_model_value"
    expected_unbound_energy = -3.12 
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_model_compact( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound_model is set to compact."""
    dpf_stem = "1pgp_unbound_model_compact"
    expected_unbound_energy =  0.00 #@FixMe 3/2009 do not know how to calc this
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_model_bound( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound_model is set to bound."""
    dpf_stem = "1pgp_unbound_model_bound"
    expected_unbound_energy = -3.23
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_set0_test( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound is set to 0 in the DPF."""
    dpf_stem = "1pgp_unbound_set0"
    expected_unbound_energy = 0.
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_set10_test( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound is set to 10 in the DPF."""
    dpf_stem = "1pgp_unbound_set10"
    expected_unbound_energy = 10.
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_1_4_50_test( AutoDock4_energy_test ):
    """Test that autodock 4.2 works include_1_4 set to 50 in the DPF."""
    dpf_stem = "1pgp_1_4_50"
    #print "in new ga_only test"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_binding_energy = -4.08
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________


class AutoDock4_1pgp_1_4_100_test( AutoDock4_energy_test ):
    """Test that autodock 4.2 works include_1_4 set to 100 in the DPF."""
    dpf_stem = "1pgp_1_4_100"
    #print "in new ga_only test"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    #expected_binding_energy = -4.08
    expected_binding_energy = -3.49
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

#print "docked[binding_energy]=", docked['binding_energy']
class AutoDock4_1pgp_ligrand_ga_only_energy_test( AutoDock4_energy_test ):
    """Test that autodock 4.2 works when ligand is randomized within autodock as set in the DPF.
                 and that expected energy is found"""
    dpf_stem = "1pgp_ligrand_ga_only"
    expected_binding_energy = 1.39
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________


class AutoDock4_1pgp_ga_only_energy_test( AutoDock4_energy_test ):
    """Test that autodock 4.2 works when ga_only is set in the DPF."""
    dpf_stem = "1pgp_ga_only"
    #print "in new ga_only test"
        #expected_binding_energy = +843.59     -5.87
        #expected_intermol_energy = +21.00  -6.17
        #expected_internal_energy = +820.50 -3.23
    expected_binding_energy = -4.08
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________


class AutoDock4_1pgp_flexres_energy_test( AutoDock4_energy_test ):
    """Test that autodock 4.2 works when flexres is set in the DPF."""
    dpf_stem = "1pgp_flexres"
    #print "in new flexres test"
        #expected_binding_energy =  -4.22
        #expected_intermol_energy = -5.93
        #expected_internal_energy = -2.18
    expected_binding_energy = -4.33
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_flexres_energy_1_4_50_test( AutoDock4_energy_test ):
    """Test that autodock 4.2 works when flexres is set in the DPF."""
    dpf_stem = "1pgp_flexres_1_4_50"
    #print "in new flexres test"
        #expected_binding_energy =  -4.22
        #expected_intermol_energy = -5.93
        #expected_internal_energy = -2.18
    #expected_binding_energy = -4.72
    expected_binding_energy = -4.33
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock4_1pgp_flexres_energy_1_4_100_test( AutoDock4_energy_test ):
    """Test that autodock 4.2 works when flexres is set in the DPF."""
    dpf_stem = "1pgp_flexres_1_4_100"
    #print "in new flexres test"
        #expected_binding_energy =  -4.22
        #expected_intermol_energy = -5.93
        #expected_internal_energy = -2.18
    expected_binding_energy = -4.33
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________


class AutoDock4_1pgp_ga_smooth0_energy_test( AutoDock4_energy_test ):
    """Test that autodock 4.2 gives expected internal energy when smooth is set to 0 in the DPF."""
    dpf_stem = "1pgp_smooth0"
    #restores former values (less negative)
    expected_binding_energy =  -5.87
    expected_intermol_energy = -6.17
    expected_internal_energy = -1.78
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________




if __name__ == '__main__':
    #  This syntax lets us run all the tests,
    #  or conveniently comment out tests we're not interested in.
    #  NOTE:  Remember to add new TestCase class names to the list "test_cases"
    test_cases = [
        # tests for .dpf extension:
        'AutoDock4_1pgp_no_extension',
        'AutoDock4_1pgp_wrong_extension',
        'AutoDock4_1pgp_two_extensions',
        # simple tests:
        'AutoDock4_1pgp_timepid_test',
	'AutoDock4_1pgp_seed01_test',
	'AutoDock4_1pgp_seed37_test',
	'AutoDock4_1pgp_seedinttime_test',
        'AutoDock4_1pgp_epdb_flexres_intelec_test',
        'AutoDock4_1pgp_epdb_flexres_nointelec_test',
        'AutoDock4_1pgp_epdb_etables_test',
        'AutoDock4_1pgp_numoutlev_test',
        'AutoDock4_1pgp_symoutlev_test',
	'AutoDock4_1pgp_no_outlev_test',
	'AutoDock4_1pgp_illegal_outlev_test',
        'AutoDock4_1pgp_illegal_symoutlev_test',
	'AutoDock4_1pgp_mixed_basic_analysis_test',
	'AutoDock4_1pgp_mixed_adt_analysis_test',
        'AutoDock4_1pgp_ligand_types_map_mismatch',
        'AutoDock4_1pgp_map_truncated',
        'AutoDock4_1pgp_illegal_keyword_test',
        'AutoDock4_1pgp_no_elecmap_test',
        'AutoDock4_1pgp_no_desolvmap_test',
        'AutoDock4_1pgp_no_elec_desolv_maps_test',
        'AutoDock4_1pgp_too_many_ligand_types_test',
        'AutoDock4_1pgp_too_many_torsions',
        'AutoDock4_1pgp_just_right_number_torsions',
        'AutoDock4_1pgp_two_ligands_test',
        'AutoDock4_1pgp_two_mapsets_test',
        'AutoDock4_1pgp_unbound_set_illegal_test',
        'AutoDock4_1pgp_unbound_model_illegal_test', #1
	'AutoDock4_1pgp_intnbpreps_test',
	'AutoDock4_1pgp_intnbpreps_toofew_test',
	'AutoDock4_1pgp_intnbpreps_illegal_atom_test',
	'AutoDock4_1pgp_intnbcoeffs_test',
        'AutoDock4_1pgp_ga_select_tournament_test',
        'AutoDock4_1pgp_ga_select_linear_ranking_test',
        #'AutoDock4_1pgp_rmsmode_heavy_atoms_only_test',
        'AutoDock4_1pgp_rmsmode_heavy_atoms_only_noH_test',
        'AutoDock4_1pgp_rmsmode_all_heavy_atom_pairs_only_test',
        'AutoDock4_1pgp_rmsmode_unique_heavy_atom_pairs_only_test',
        'AutoDock4_1pgp_simanneal_small', 
        'AutoDock4_1pgp_simanneal_small_geometric_schedule', 
        #next dpf sets tran0,quaternion0,dihe0 to random
        'AutoDock4_1pgp_ligrand_ga_only_test', 
        'AutoDock4_1pgp_ga_only_test',
	'AutoDock4_1pgp_output_pop_file_test',
        # tests using flexible residues 
        'AutoDock4_1pgp_flexres_test',
        ## tests which check for specific value
        'AutoDock4_1pgp_test',
        'AutoDock4_1pgp_smaller_test',
## next test not installed yet as I do not have a way I know of to pass
## different expected values to the superclass - M Pique  July 2012
	# 'AutoDock4_1pgp_intelec_off_test',
        'AutoDock4_1pgp_no_parameter_file_test',
	'AutoDock4_1pgp_gals_set_sw1_test',
	'AutoDock4_1pgp_no_set_psw1_test',
	'AutoDock4_1pgp_no_set_ga_test',
	'AutoDock4_1pgp_gals_use_defaults_test',
	'AutoDock4_1pgp_simanneal_10runs_test',
	'AutoDock4_1pgp_lsonly_10runs_test',
	'AutoDock4_1pgp_gals_10runs_test',
        #'AutoDock4_1pgp_ga_only_value_test',
        ## tests for unbound values 
        'AutoDock4_1pgp_unbound_default_test',
        'AutoDock4_1pgp_unbound_set0_test',
        'AutoDock4_1pgp_unbound_set10_test',
        # tests for unbound_model choices
        'AutoDock4_1pgp_unbound_model_compact',
        'AutoDock4_1pgp_unbound_model_bound',
        'AutoDock4_1pgp_unbound_model_extended',
        'AutoDock4_1pgp_unbound_compute_unbound_extended',
        'AutoDock4_1pgp_unbound_model_value',
        # tests for ga_only and energy
        'AutoDock4_1pgp_ligrand_ga_only_energy_test',
        'AutoDock4_1pgp_ga_only_energy_test',
        'AutoDock4_1pgp_ga_smooth0_energy_test',
        # tests for energy with flexible residues 
        'AutoDock4_1pgp_flexres_energy_test',
        'AutoDock4_1pgp_flexres_energy_1_4_50_test',
        'AutoDock4_1pgp_flexres_energy_1_4_100_test',
        # tests for energy with 1_4_interactions
        'AutoDock4_1pgp_1_4_50_test',
        'AutoDock4_1pgp_1_4_100_test',
        # tests for setting tran0 from about
        'AutoDock4_1pgp_about_only_test', 
        'AutoDock4_1pgp_multipleabouts_test', 
        # tests for setting 'about' automatically
        'AutoDock4_1pgp_about_auto_simanneal', 
        'AutoDock4_1pgp_about_auto_ga_only', 
        'AutoDock4_1pgp_about_auto_local_only', 
        'AutoDock4_1pgp_about_auto_gals', 
	# tests that take a long time to run 
	'AutoDock4_1pgp_ga_run_maxruns_test',
	'AutoDock4_1pgp_ga_run_maxruns_analysis_test',
	'AutoDock4_1pgp_mixed_maxruns_test',
	'AutoDock4_1pgp_mixed_maxruns_analysis_test',
	'AutoDock4_1pgp_mixed_maxruns_basic_analysis_test',
	'AutoDock4_1pgp_mixed_maxruns_adt_analysis_test',
	'AutoDock4_1pgp_overmaxruns_test',
    ]
    unittest.main( argv=( [__name__ ,] + test_cases ) )
    #  The call "unittest.main()" automatically runs all the TestCase classes in
    #  alphabetical order; calling with argv=([]), lets us specify the order.
    #  NOTE: "unittest.main()" saves us having to remember to add new tests to the 
    #  list of test cases.
    #unittest.main()
    #  For verbose output, use this:
    #unittest.main( argv=( [__name__, '-v'] + test_cases ) )
