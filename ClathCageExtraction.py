###############################################################################
# Program to Extract Cages from the clathrate hydrate unit cell
# Written by Joshua Cantin on 30 Apr 2013
###############################################################################
#Need:
#	-xyz file with locations of oxygen and hydrogen atoms for the clathrate hydrate unit cell
# 	-xyz file with the cage centre locations
#       -Clath_xyz.py, a module upon which this code depends
###############################################################################
# To run: pydev ClathCageExtraction.py [xyz_atomFile] [xyz_cageFile] [structure]
###########
# [xyz_atomFile]: Contains the xyz coordinates of O and H atoms
# [xyz_cageFile]: Contains the xyz coordinates of the centres of each cage
# [structure]:    Indicates the structure of the clathrate hydrate
#                     Can be:
#                        sII - for the sII structure
#                        sI  - for the sI structure
###########
# Output: Outputs one xyz file for each cage, the supercell (3X3) xyz file,
#         the orginal reference files and a test data file in a subdirectory
#         labelled with the cell structure and the date and time.
###############################################################################################
# NOTE: The code makes the following assumptions:
# 1. All of the unit cells have their origins as the back left corner of the cell, as seen in
#    Takeuchi, et al. Water proton configurations in structures I, II, and H clathrate hydrate unit cells. J. Chem. Phys. 2013, 138, 124504.
# 2. The cages are empty and the closest oxygen atoms to the cage centre are those that directly form the cage in question.

#References:
#     Takeuchi, et al. Water proton configurations in structures I, II, and H clathrate hydrate unit cells. J. Chem. Phys. 2013, 138, 124504.
#     Ratcliffe, C. I.; Ripmeester, J. A. Proton and carbon-13 NMR studies on carbon dioxide hydrate. J. Phys. Chem. 1986, 90, 1259-1263.

import numpy as np
import Scientific.Geometry as Geo
from Scientific.IO.TextFile import TextFile
import string
import sys
import os
import datetime
import Clath_xyz

###############################################################################
#Function Definitions
###############################################################################

def getCageAtomVecs(oxygenAtoms, hydrogenAtoms, OH1_correlations, OH2_correlations, Cage_Os, cageCenters):

    ###############################################################################
    # Get cage vectors given which oxygens belong to the cage.
    ###############################################################################

    Cage_O_vecs = [[]]
    Cage_H_vecs = [[]]
    Cage_OH1_cor = [[]]
    Cage_OH2_cor = [[]]

    # Get vectors referenced back to the origin
    for cageNo in range(0, len(Cage_Os)):
        i = 0
        for oxyNo in Cage_Os[cageNo]:
            Cage_O_vecs[cageNo].append(oxygenAtoms[oxyNo] - cageCenters[cageNo])
            Cage_H_vecs[cageNo].append(hydrogenAtoms[OH1_correlations[oxyNo]] - cageCenters[cageNo])
            Cage_OH1_cor[cageNo].append(i)
            Cage_H_vecs[cageNo].append(hydrogenAtoms[OH2_correlations[oxyNo]] - cageCenters[cageNo])
            Cage_OH2_cor[cageNo].append(i+1)
            i += 2
            
        Cage_O_vecs.append([])
        Cage_H_vecs.append([])
        Cage_OH1_cor.append([])
        Cage_OH2_cor.append([])

    return {"Cage_O_vecs" : Cage_O_vecs, "Cage_H_vecs" : Cage_H_vecs, "Cage_OH1_cor" : Cage_OH1_cor, "Cage_OH2_cor" : Cage_OH2_cor}

###############################################################################
#Start of Program
###############################################################################


lenOH = 1.1 #Angstroms; assumed maximum length of the OH bond

#Average cage radii, taken from Sloan, E. D. Clathrate hydrates of natural gases; Chemical industries; 3rd ed.; CRC Press: Boca Raton, FL, 2008. 
sI_cageRadii =  {"Ar(Small)": 3.95, "Kr(Large)": 4.33} #Angstroms
sII_cageRadii = {"Ar(Small)": 3.91, "Kr(Large)": 4.73} #Angstroms
cageRadii = {"sI": sI_cageRadii, "sII": sII_cageRadii}

#The number of waters in each cage, taken from Sloan, E. D. Clathrate hydrates of natural gases; Chemical industries; 3rd ed.; CRC Press: Boca Raton, FL, 2008.
sI_cageWaterCount = {"Ar(Small)": 20, "Kr(Large)": 24}
sII_cageWaterCount = {"Ar(Small)": 20, "Kr(Large)": 28}
cageWaterCount = {"sI": sI_cageWaterCount, "sII": sII_cageWaterCount}

#The Lattice Constants for the sI and sII unit cells:
# According to Takeuchi, et al.:
#  Structure      Lattice Constant (A)
#    sII               17.31 (Takeuchi, et al.)
#    sI                12.03 (Takeuchi, et al. say that Ratcliffe and Ripmeester have this, but I didn't find it after a quick look)
latticeConst = {"sI": 12.03, "sII": 17.31}

xyz_atomFilename = sys.argv[1]
xyz_cageFilename = sys.argv[2]
cellStructure = sys.argv[3]
cellSize = latticeConst[cellStructure]
#cellSize = float(sys.argv[4])

###############################################################################
#Read in the O and H atom locations
###############################################################################
out_dict = Clath_xyz.readWaterAtoms(xyz_atomFilename)

Hyd = out_dict["H"]
Oxy = out_dict["O"]

print "Read in atom positions. There are %d H atoms and %s O atoms." % (len(Hyd), len(Oxy))


###############################################################################
#Read in the Cage Centre locations
###############################################################################
large_cages = []
small_cages = []

#Open input file
xyz_file = TextFile(xyz_cageFilename)

#Ignore first two lines of the file
junk = xyz_file.readline()
junk = xyz_file.readline()

while 1:
    line = xyz_file.readline()
    if not line: break

    x, y, z = map(float, string.split(line)[1:]) #ignore atom label
    cage_type = string.split(line)[0] #get the atom type

    #Convert coordinates based on atom type

    if (cage_type == "Ar(Small)"):

        #Convert coordinates to vectors for the Small Cages
        small_cages.append(Geo.Vector(x,y,z))

    elif (cage_type == "Kr(Large)"):

        #Convert coordinates to vectors for the Large Cages
        large_cages.append(Geo.Vector(x,y,z))

    else:

        raise CorrelationError("Can't parse file. Encountered strange cage centre type: " + cage_type + " Only 'Ar(Small)' and 'Kr(Large)' accepted.")

xyz_file.close()

print "Read in Cage Centre locations. There are %d small and %d large cages." % (len(small_cages), len(large_cages))


###############################################################################
#Generate 27 duplicates of the unit cell to make a 3X3 cube of unit cells
#     - Do this so that all of the cages in the central unit cell are fully formed 
###############################################################################
Hyd_cell = []
Oxy_cell = []

for i in range(-1, 2):
    for j in range(-1, 2):
        for k in range(-1, 2):
            #print "###################################################################################"
            #print (i,j,k)
            for l in range(0, len(Hyd)):
                Hyd_cell.append(Hyd[l] + (Geo.Vector(i,j,k) * cellSize))
                #print Hyd[l] + (Geo.Vector(i,j,k) * cellSize)
            for m in range(0, len(Oxy)):
                Oxy_cell.append(Oxy[m] + (Geo.Vector(i,j,k) * cellSize))
                #print Oxy[m] + (Geo.Vector(i,j,k) * cellSize)

print "Made unit cell duplicates."


###############################################################################
#Determine which H's go with which O's
###############################################################################
out_dict = Clath_xyz.correlateWaters(Oxy_cell, Hyd_cell, lenOH)

OH1_all_cor = out_dict["OH1_all_cor"]
OH2_all_cor = out_dict["OH2_all_cor"]

print "Correlated O's and H's." #This piece takes the longest amount of time so far.


###############################################################################
#Calculate internal water molecule vectors OH1, OH2, and HH
###############################################################################
OH1 = []
OH2 = []
HH = []
for waterNo in range(0,len(Oxy_cell)):
    OH1.append(Hyd_cell[OH1_all_cor[waterNo]] - Oxy_cell[waterNo])
    OH2.append(Hyd_cell[OH2_all_cor[waterNo]] - Oxy_cell[waterNo])
    HH.append(Hyd_cell[OH2_all_cor[waterNo]] - Hyd_cell[OH1_all_cor[waterNo]])

print "Calculated bond vectors."


###############################################################################
#Determine the O atoms in each cage
###############################################################################

#Determine the O atoms by taking the correct number of O atoms nearest to the cage centre.
#Double check by calculating the mean radius and comparing it to lit values (Sloan and Koh)
#and by calcuating the average O atom position and comparing it to the cage centre.

##############
# Small Cages
##############

out_dict = Clath_xyz.extractCages(Oxy_cell, small_cages, cellStructure, "small", cageWaterCount, cageRadii)

smallCage_Os = out_dict["Cage_Os"]
smallCage_TestsStr = out_dict["CageTestsStr"]
smallCageCentres_vecAvg = out_dict["cageCentreFromAvg"]

for testStr in smallCage_TestsStr:
    print testStr

out_dict = getCageAtomVecs(Oxy_cell, Hyd_cell, OH1_all_cor, OH2_all_cor, smallCage_Os, smallCageCentres_vecAvg)

smallCage_O_vecs = out_dict["Cage_O_vecs"]
smallCage_H_vecs = out_dict["Cage_H_vecs"]
smallCage_OH1_cor = out_dict["Cage_OH1_cor"]
smallCage_OH2_cor = out_dict["Cage_OH2_cor"]
    
##############
# Large Cages
##############

out_dict = Clath_xyz.extractCages(Oxy_cell, large_cages, cellStructure, "large", cageWaterCount, cageRadii)

largeCage_Os = out_dict["Cage_Os"]
largeCage_TestsStr = out_dict["CageTestsStr"]
largeCageCentres_vecAvg = out_dict["cageCentreFromAvg"]

for testStr in largeCage_TestsStr:
    print testStr

out_dict = getCageAtomVecs(Oxy_cell, Hyd_cell, OH1_all_cor, OH2_all_cor, largeCage_Os, largeCageCentres_vecAvg)

largeCage_O_vecs = out_dict["Cage_O_vecs"]
largeCage_H_vecs = out_dict["Cage_H_vecs"]
largeCage_OH1_cor = out_dict["Cage_OH1_cor"]
largeCage_OH2_cor = out_dict["Cage_OH2_cor"]
    
print "Cages extracted."

#print largeCage_Os

#print large_cages

###############################################################################
#Print the Cages to .xyz files in a subdirectory labelled with structure and date/time
###############################################################################

#Get the date
date_time = datetime.datetime.now()

#Make the subdirectory
subdirectoryName = cellStructure + "_cages_" + date_time.strftime("%Y-%m-%d_%H%M")
os.system("mkdir %s" % subdirectoryName)

# Output small cage xyz files in the subdirectory
for cageNo in range(0, len(small_cages)):

    cageFilename = subdirectoryName + "/" + cellStructure + "_smallCage" + str(cageNo) + ".xyz"
    cageHeader = cellStructure + " Clathrate Hydrate " + "Small" + " Cage " + str(cageNo) +"; x y z coordinates in Angstroms"

    Clath_xyz.saveWatersToXYZ(smallCage_O_vecs[cageNo], smallCage_H_vecs[cageNo], smallCage_OH1_cor[cageNo], smallCage_OH2_cor[cageNo], cageFilename, cageHeader)

# Output large cage xyz files in the subdirectory
for cageNo in range(0, len(large_cages)):

    cageFilename = subdirectoryName + "/" + cellStructure + "_largeCage" + str(cageNo) + ".xyz"
    cageHeader = cellStructure + " Clathrate Hydrate " + "Large" + " Cage " + str(cageNo) +"; x y z coordinates in Angstroms"

    Clath_xyz.saveWatersToXYZ(largeCage_O_vecs[cageNo], largeCage_H_vecs[cageNo], largeCage_OH1_cor[cageNo], largeCage_OH2_cor[cageNo], cageFilename, cageHeader)

print "Cage .xyz files output to ./" + subdirectoryName + "/"

###############################################################################
#Save Cage Test Data
###############################################################################
cageTestFilename = subdirectoryName + "/" + cellStructure + "_cageData.txt"

cageTestFile = open(cageTestFilename, 'w')

cageTestFile.write("Data concerning the tests performed while extracting the cages from the %s unit cell. \nOriginal Unit Cell from file %s and original cage centre data from file %s."
                   % (cellStructure, xyz_atomFilename, xyz_cageFilename) 
                   + "\n")

cageTestFile.write("##############################################################################################################" + "\n")

for testStr in smallCage_TestsStr:
    cageTestFile.write(testStr + "\n")
    
cageTestFile.write("##############################################################################################################" + "\n")

for testStr in largeCage_TestsStr:
    cageTestFile.write(testStr + "\n")

print "Cage Test data written to ./" + cageTestFilename


###############################################################################
#Print the Whole SuperCell to an .xyz file for input into VMD
###############################################################################

superCellFilename = subdirectoryName + "/" + cellStructure + "_superCell.xyz"
superCellHeader = "3X3 unit cells; x y z coordinates in Angstroms"

Clath_xyz.saveWatersToXYZ(Oxy_cell, Hyd_cell, OH1_all_cor, OH2_all_cor, superCellFilename, superCellHeader)

print ("XYZ Supercell data written to ./" + superCellFilename)

###############################################################################
#Copy over original files for reference purposes
###############################################################################
origin = xyz_atomFilename
destination = "./" + subdirectoryName + "/" + origin
os.system("cp " + origin + " " + destination)

origin = xyz_cageFilename
destination = "./" + subdirectoryName + "/" + origin
os.system("cp " + origin + " " + destination)

print "%s and %s copied to ./%s for reference purposes." % (xyz_atomFilename, xyz_cageFilename, subdirectoryName)
