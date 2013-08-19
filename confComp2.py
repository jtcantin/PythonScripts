###############################################################################################
# Configuration Comparer - confComp2.py
###############################################################################################
# This program reads in the files and then compares the two configurations
###############################################################################################
# To run: pydev confComp.py [xyz_atomFile] [xyz_CMEAFile]
###########
# [xyz_atomFile]: Contains the xyz coordinates of each atom of each molecule
# [xyz_CMEAFile]: Contains the xyz coordinates of the centre of mass (CM) and Euler Angles (EA) of each molecule
###########
# Output: Returns the RMS difference between the two different configurations
###############################################################################################

from MMTK import *
import sys
import ClathIO as CIO
from Scientific.IO.TextFile import TextFile
import numpy as np

universe = InfiniteUniverse()
#universe = OrthorhombicPeriodicUniverse((2,2,2))

xyz_atomFile = TextFile(sys.argv[1])
xyz_CMEAFile = TextFile(sys.argv[2])

#Set positions of O and H atoms from original xyz file
CIO.readClathrate_xyz(xyz_atomFile, universe)

#universe.view()

conf_xyz = copy(universe.configuration())

# Empty the universe
#universe = OrthorhombicPeriodicUniverse((2,2,2))
#universe = InfiniteUniverse()

#Set water molecule positions from ZY'Z'' Euler Angles and Centre of Mass data
CIO.readClathrate104_ZYZ_oneBead(xyz_CMEAFile, universe)

#Set water molecule positions from zxz Euler Angles and Centre of Mass data
#CIO.readClathrate104(xyz_CMEAFile, universe)

#universe.view()

conf_ZYZ = copy(universe.configuration())

# Determine the RMS difference between the locations of each atom
# NOTE: The below code assumes that the .xyz file is written as O H1 H2 for each water molecule
sumdiffsq = 0
maxdiffsq = 0

first_water_index = next(i for i in range(0,len(universe.objectList())) if universe[i].fullName()=="water")


for objNo in range(0, len(universe.objectList())):
    if objNo < (first_water_index):
         if (universe[objNo].fullName() == "oxygen"):
             diff = universe[objNo].position() - universe[first_water_index+int(np.floor(objNo/3))].O.position()
             diffsq = diff*diff
             sumdiffsq += diffsq
             maxdiffsq = max(maxdiffsq, diffsq)
                          
         elif (universe[objNo].fullName() == "hydrogen" and universe[objNo-1].fullName() == "hydrogen"):
             diff = universe[objNo].position() - universe[first_water_index+int(np.floor(objNo/3))].H2.position()
             diffsq = diff*diff
             sumdiffsq += diffsq
             maxdiffsq = max(maxdiffsq, diffsq)
                          
         elif (universe[objNo].fullName() == "hydrogen" and universe[objNo+1].fullName() == "hydrogen"):
             diff = universe[objNo].position() - universe[first_water_index+int(np.floor(objNo/3))].H1.position()
             diffsq = diff*diff
             sumdiffsq += diffsq
             maxdiffsq = max(maxdiffsq, diffsq)
                        
         else: print "Error"

    pass
print first_water_index
meansumdiffsq = sumdiffsq/(first_water_index) # The first_water_index is the same as the number of H and O atoms combined, since the list starts at 0 and not 1
rootmeansumdiffsq = np.sqrt(meansumdiffsq)
print "RMS Diff = %.5e nm" % (rootmeansumdiffsq)
print "Max |Diff| = %.5e nm" % (np.sqrt(maxdiffsq))



