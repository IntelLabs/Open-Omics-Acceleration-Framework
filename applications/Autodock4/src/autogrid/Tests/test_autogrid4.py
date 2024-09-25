#! /usr/bin/env python
# 
#
# $Id: test_autogrid4.py,v 1.18 2013/11/14 22:16:02 mp Exp $
#
"""
Test AutoGrid.
"""

#______________________________________________________________________________
#
import sys
import os
import types
import unittest
import getopt
import subprocess
from string import split, strip

#______________________________________________________________________________
#
# Global variables

autogrid_executable = "../autogrid4" # where the AutoGrid executable resides
gpf_directory = '.' # where the input GPF files reside
test_output_directory = '.' # where the GLG files will be written

try:
    opts, argv = getopt.getopt(sys.argv[1:], "d:g:e:o:",
    ["gpf-directory=","gpf-directory=","executable=","test-output-directory="]) 
except getopt.GetoptError, v:
    usage() 
    sys.exit(2)

for o,a in opts:
    if o in ("-d", "--gpf-directory"):
        gpf_directory = a
    if o in ("-g", "--gpf-directory"):
        gpf_directory = a
    if o in ("-e", "--executable"):
        autogrid_executable = a
    if o in ("-o","--test-output-directory"):
        test_output_directory = a


built_maps = False
#no parameter library keyword 
built_maps_no_parameter_library = False
#parse all the receptor types from the receptor input file
built_maps_no_receptor_types = False 
#parse two extra receptor types from the receptor input file
built_maps_minus_two_types = False  
#no receptor types, ligand_types keyword preceeds receptor filename
built_maps_ligand_types_before_receptor = False 
#missing_map for ligand_types keyword
built_missing_map = False

#______________________________________________________________________________

def usage():
    """Print out the usage of this command."""
    print """Usage:  python test_autogrid4.py [-d <string>] [-e <string>] [-o <string>]

where:
    -d, --gpf-directory
        specifies the directory containing the GPFs to be tested;
        this flag is optional; default is '.'
    -e, --executable
        specifies the path to the AutoGrid executable to be tested;
        this flag is optional; default is '../autogrid4'
    -o, --test-output-directory
        specifies the directory where the output GLGs will be written;
        this flag is optional; default is '.'

NOTE:  these may be relative to the directory where this script was invoked.
"""

#______________________________________________________________________________

def run_AutoGrid( gpf_filename, glg_filename ):
    """Launch AutoGrid, using the specified AutoGrid executable and GPF,
    create the specified GLG, and trap all the outputs from standard output
    and standard error."""
    gpf = gpf_directory + os.sep + gpf_filename
    glg = test_output_directory + os.sep + glg_filename
    rm(glg)
    command = [ autogrid_executable, '-p', gpf, '-l', glg ]
    print '\nRunning ' + autogrid_executable + ' using GPF "'+gpf+'", saving results in "'+glg+'":'
    try:
	rc = subprocess.call( command )
	#print 'autogrid returned ', rc  # DEBUG
        return find_success_in_GLG(glg_filename)
    except OSError, e:
        print "\nUnable to run " + autogrid_executable + " :", e
        return False

#______________________________________________________________________________

def find_success_in_GLG( glg_filename ):
    """Open the AutoGrid GLG, and look for the string "Successful Completion"
    in the last 10 lines of the file."""
    glg = test_output_directory + os.sep + glg_filename
    try:
        fptr = open( glg )
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
    """Remove single file"""
    # no regular expression wildcards, e.g. * or ?, allowed
    if os.access(filename, os.F_OK) :
       try:
          os.remove( filename )
       except OSError,e:
         print 'unable to remove '+filename

#______________________________________________________________________________

class AutoGrid_hsg1_sm_test(unittest.TestCase):
    
    def setUp(self):
        """Set up for autogrid4 tests.
        Locate the autogrid binary now during setUp.
        """
        global built_maps
        self.autogrid = autogrid_executable

        if not built_maps:
            # Make sure you remove all the products of AutoGrid from any previous tests.
            #rm("hsg1_sm.*map*")
            built_maps = run_AutoGrid('hsg1_sm.gpf', 'hsg1_sm.glg')

    def tearDown(self):
        #??remove map here
        pass

   # def run_cmd(self, cmd_str):
   #   global built_maps
   #     (i,o,e) = os.popen3(cmd_str) # trap all the outputs
   #     #print 'calculating maps...'
   #     os.wait() # for the child process to finish
   #     built_maps = True

    def compare_autogrid4_maps(self, stem, maptype, precision=2):
        """Compare autogrid4_maps from C vs python
            allow .225 vs .226 to pass using 'precision'
            allow .526 vs .524 to pass using 0.5 'percentage' check
        """
        # read in 'c' and 'py' maps
        c_autogrid4_map = test_output_directory + os.sep + stem + '.' + maptype +'.map'
        fptr = open(c_autogrid4_map)
        c_autogrid4_lines = fptr.readlines()
        fptr.close()
        #read in python AutoGrid4 results
        py_autogrid4_map = gpf_directory + os.sep + 'py_' + stem +'.'+ maptype+'.map'
        fptr = open(py_autogrid4_map)
        py_autogrid4_lines = fptr.readlines()
        fptr.close()
        # skip starting lines
        self.assertEqual(len(c_autogrid4_lines), len(py_autogrid4_lines))
        for atomL, py_atomL in zip(c_autogrid4_lines[6:], py_autogrid4_lines[6:]):
            c_num = float(strip(atomL))
            py_num = float(strip(py_atomL))
            if abs(c_num)>.30:
                # use percentage difference for large(r) values:
                #  allowing up to 0.5% error 
                #if ((abs(c_num-py_num)/(c_num))*100.<.5)==False:
                #    print  "c_num=", c_num," py_num=", py_num, "test=", (abs(c_num-py_num)/(c_num))*100.
                self.assertEquals(((abs(c_num-py_num)/(c_num))*100.)<.5, True)
            else:
                #but use precision only for smaller values
                cutoff = 10**(-precision)
                self.assertEqual(abs(c_num-py_num)<cutoff, True)
                #self.assertAlmostEquals(c_num, py_num, precision)
        # clean up
        #rm(c_autogrid4_map)

    def test_hsg1_estat(self):
        """test_hsg1 estat case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'e')

    def test_hsg1_dsolv(self):
        """test_hsg1 dsolv case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'd')

    def test_hsg1_A(self):
        """test_hsg1 A case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'A')

    def test_hsg1_C(self):
        """test_hsg1 C case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'C')
 
    def test_hsg1_N(self):
        """test_hsg1 N case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'N')

    def test_hsg1_HD(self):
        """test_hsg1 HD case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'HD')

    def test_hsg1_NA(self):
        """test_hsg1 NA case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'NA')

    def test_hsg1_OA(self):
        """test_hsg1 OA case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'OA')


class AutoGrid_hsg1_sm_no_parameter_library_test(AutoGrid_hsg1_sm_test):

    def setUp(self):
        """Set up for autogrid4 tests.  Locate the autogrid binary now during setUp."""
        global built_maps_no_parameter_library
        self.autogrid = autogrid_executable
        if not built_maps_no_parameter_library:
            # Make sure you remove all the products of AutoGrid from any previous tests.
            #rm("hsg1_sm.*map*")
            built_maps_no_parameter_library = run_AutoGrid( 'hsg1_sm_no_parameter_file.gpf', 'hsg1_sm.glg' )


class AutoGrid_hsg1_sm_no_receptor_types_test(AutoGrid_hsg1_sm_test):

    def setUp(self):
        """Set up for autogrid4 tests.  Locate the autogrid binary now during setUp."""
        global built_maps_no_receptor_types
        self.autogrid = autogrid_executable
        if not built_maps_no_receptor_types:
            # Make sure you remove all the products of AutoGrid from any previous tests.
            #rm("hsg1_sm.*map*")
            built_maps_no_receptor_types = run_AutoGrid( 'hsg1_no_receptor_types.gpf', 'hsg1_sm.glg' )


class AutoGrid_hsg1_sm_minus_two_types_test(AutoGrid_hsg1_sm_test):

    def setUp(self):
        """Set up for autogrid4 tests.  Locate the autogrid binary now during setUp."""
        global built_maps_minus_two_types
        self.autogrid = autogrid_executable
        if not built_maps_minus_two_types:
            # Make sure you remove all the products of AutoGrid from any previous tests.
            #rm("hsg1_sm.*map*")
            built_maps_minus_two_types = run_AutoGrid( 'hsg1_no_receptor_types.gpf', 'hsg1_sm.glg' )


class AutoGrid_ligand_types_before_receptor_test(AutoGrid_hsg1_sm_test):

    def setUp(self):
        """Set up for autogrid4 tests. Locate the autogrid binary now during setUp."""
        global built_maps_ligand_types_before_receptor
        self.autogrid = autogrid_executable
        if not built_maps_ligand_types_before_receptor:
            # Make sure you remove all the products of AutoGrid from any previous tests.
            #rm("hsg1_sm.*map*")
            built_maps_ligand_types_before_receptor = run_AutoGrid( 'hsg1_ligand_types_before_receptor.gpf', 'hsg1_sm.glg')

#------------------------------------------------------------------
# tests contributed by Stefano Forli
#------------------------------------------------------------------
class AutoGrid_simple_test(unittest.TestCase):
    """ Base Class for AutoGrid testing."""
    gpf_stem = "BaseClass"
    computed = False

    def setUp(self):
        """ Set up for autogrid4 tests. Locate the autogrid binary now during setUp."""
        self.glg_filename = "test_" + self.gpf_stem + ".glg"
        self.computed = run_AutoGrid(self.gpf_stem + ".gpf", self.glg_filename)
        print "after call to run_AutoGrid"

    def test_glg_exists(self):
        #Check that autogrid calculation finished and a new GLG has been computed
        if (self.expected_outcome == True):
            print "Testing that GLG exists and AutoGrid successfully completed."
        else:
            print "Testing that GLG exists and AutoGrid did not complete."
        self.assertEqual(self.computed, self.expected_outcome)

#------------------------------------------------------------------
class AutoGrid_control_test(AutoGrid_simple_test):
    gpf_stem = 'x1hpv_CASE0'  #not missing A,elecmap or dsolvmap
    expected_outcome = True

#------------------------------------------------------------------
class AutoGrid_missing_Amap_test(AutoGrid_simple_test):
    gpf_stem = 'x1hpv_CASE1'  #missing x1hpv.A.map
    expected_outcome = False

#------------------------------------------------------------------
class AutoGrid_missing_elecmap_test(AutoGrid_simple_test):
    gpf_stem = 'x1hpv_CASE2'  #missing elecmap
    expected_outcome = False

#------------------------------------------------------------------
class AutoGrid_missing_dsolvmap_test(AutoGrid_simple_test):
    gpf_stem = 'x1hpv_CASE3'  #missing dsolvmap
    expected_outcome = False

#------------------------------------------------------------------
class AutoGrid_nbpreps_test(AutoGrid_simple_test):
    gpf_stem = 'hsg1_sm_nbpreps'  # nbp_r_eps
    expected_outcome = True

#------------------------------------------------------------------
class AutoGrid_nbpreps_illegal_atom_test(AutoGrid_simple_test):
    gpf_stem = 'hsg1_sm_nbpreps_illegal_atom'  # nbp_r_eps
    # illegal atoms in nbb_r_eps statements are not considered an error
    expected_outcome = True

#------------------------------------------------------------------
class AutoGrid_nbpreps_toofew_test(AutoGrid_simple_test):
    gpf_stem = 'hsg1_sm_nbpreps_toofew'  # nbp_r_eps
    expected_outcome = False

#------------------------------------------------------------------
class AutoGrid_nbpcoeffs_test(AutoGrid_simple_test):
    gpf_stem = 'hsg1_sm_nbpcoeffs'  # nbp_coeffs
    expected_outcome = True



if __name__ == '__main__':
    test_cases = [
        'AutoGrid_hsg1_sm_test',
        'AutoGrid_hsg1_sm_no_parameter_library_test',
        'AutoGrid_hsg1_sm_no_receptor_types_test',
        'AutoGrid_hsg1_sm_minus_two_types_test',
        'AutoGrid_missing_Amap_test',
        'AutoGrid_missing_elecmap_test',
        'AutoGrid_missing_dsolvmap_test',
        'AutoGrid_control_test',
        #'AutoGrid_ligand_types_before_receptor_test', 
        #3/18/09 not sure of status of ligand_types_before_receptor test, rh
	'AutoGrid_nbpreps_test',
	'AutoGrid_nbpreps_illegal_atom_test',
	'AutoGrid_nbpreps_toofew_test',
	'AutoGrid_nbpcoeffs_test',
    ]
    unittest.main( argv=([__name__,] + test_cases))  # non-verbose output
    # optional:  for verbose output, use this:
    # unittest.main( argv=([__name__, '-v'] + test_cases)) # verbose output

# EOF
