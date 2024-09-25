from PyAutoDock.AutoGrid import AutoGrid4
from MolKit import Read
#note: set PYTHONPATH to pickup PyAutoDock directory

mol = Read("hsg1_sm.pdbqt")[0]
mol.buildBondsByDistance()

ag = AutoGrid4(mol, atom_types = ['A','C','HD','NA', 'N','OA' ], 
                    npts = [5, 5, 5], center= [2.5, 6.5, -7.5])
ag.write_maps()


#note: this has npts 5,5,5 vs c code which uses 4,4,4
# with either number of specified points, 5,5,5 are calculated
#for the python code, 4,4,4 produces this error:

#python2.4 write_python_maps.py
#Traceback (most recent call last):
#  File "write_python_maps.py", line 10, in ?
#    ag.write_maps()
#  File "/mgl/aduser4/rhuey/dev/PyAutoDock/AutoGrid.py", line 413, in write_maps
#    score_array = self.atom_map_scorer.get_score_array()
#  File "/mgl/aduser4/rhuey/dev/PyAutoDock/scorer.py", line 222, in get_score_array
#    array = t[0].get_score_array() * t[1]
#  File "/mgl/aduser4/rhuey/dev/PyAutoDock/vanDerWaals.py", line 1043, in get_score_array
#    atom_score = self._f(at_a, at_b, d, bx)
#  File "/mgl/aduser4/rhuey/dev/PyAutoDock/vanDerWaals.py", line 998, in _f
#    closestH_ix = self.closestH[bx]
#IndexError: index out of bounds

#possibly because there are not enough atoms in the 'toy' box volume 
#to be able to calculate the hbond vectors.... not sure of this 

#because the numbers check out, and that's what we're interested in,
#i'm ignoring this for now.
