###############################################################################
#Program to calculate various geometric properties of a clathrate
###############################################################################
#Program created on 15 Jan 2013 by Joshua Cantin and based on "CalcBondAngel_xyzFile_JTCv4.f"
#The OH bond lengths code is based on Matt Schmidt's ogiginal Fortran code
# Version 1 completed on 18 Jan 2013 by Joshua Cantin
# Version 2 started on 22 Jan 2013 by Joshua Cantin
#           -> to incoporate calculation of the centre of mass positions
#     Updated on 23 Jan 2013 to incorporate Euler Angles (note, still need to write the output
#           subroutine)
#     Updated on 24 Jan 2013 to output CM and Euler angle data to a .csv file
#     Changed on 09 Apr 2013 to be independent of HHO or OHH order in file.
#           Dangling bond count set as an option.
###############################################################################
###############################################################################
# Options
###############################################################################
# 1 = Normal for a full unit cell
# 2 = Only a single cage, perform:
#          Dangling bond count
###############################################################################

###############################################################################
#Assumptions that this code makes:
#1. That the OH bond is less than 1A and that no H's are closer to the O
#	then those that are part of the water molecule
#	-> This affects the determination of which H's are part of the same
#		water molecule as the O-atom; if violated, the distances and
#		other properties calculated will be invalid
#2. That there are only 3 nearest O neighbours to each O molecule
#	-> This affects the determination of the nearest nieghbours and
#          thus the angles of the edges of the clathrate (O-O-O angles);
#          if violated, these angles will be wrong
#3. That a dangling H-bond has |cos(theta)| > 0.5, where theta is the OH bond to
#       O position vector angle (approximately radial), and no non-dangling
#       H-bond has |cos(theta)| > 0.5.
#	-> This affects the COUNT of the number of dangling bonds; if violated
#          this will cause an incorrect counting of the number of dangling
#          H-bonds
#4. That the OH bond pointing toward the O of another water molecule is no more
#       than 1 deg off of being linear (ie. the H-O-O angle is less than 1 deg)
#       -> This will affect the determination of WHERE each OH bond points
#          (ie. to which O or if it is dangling); if violated, the OH bond
#          direction will be wrong as well as whether it is dangling or not.
###############################################################################


import numpy as np
import Scientific as Sci
import Scientific.Geometry as Geo
import sys
import string
import os

###############################################################################
#Custom Error Classes
###############################################################################

class CorrelationError(Exception): pass
        

###############################################################################
#Start of Program
###############################################################################

#Hydrogen and Oxygen masses
#Data Taken from: Standard Atomic Weights (2009). In CRC Handbook of Chemistry
#                 and Physics, 93rd ed.; CRC Press, 2012-2013.
#Atomic Weights (I have two different masses for H in case H-O-D is studied)
M_HYD1 = 1.008 #g/mol
M_HYD2 = 1.008 #g/mol
M_OXY = 15.999 #g/mol

N_a = 6.022141E23 #/mol Taken from:CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL
                  #     PHYSICAL CONSTANTS: 2010. In CRC Handbook of Chemistry
                  #     and Physics, 93rd ed.; CRC Press, 2012-2013.

m_HYD1 = M_HYD1/N_a/1000 #kg
m_HYD2 = M_HYD2/N_a/1000 #kg
m_OXY = M_OXY/N_a/1000 #kg

print "Program Executing"

input_file = sys.argv[1] 
output_file = sys.argv[2]

#input_file = "51262cageTest.xyz"
#output_file = "51262Geo.csv"

#Strip output file extension and then replace it with _CM.xyz
CM_output_file = output_file[:(string.rfind(output_file, "."))] + "_CM.xyz"
CM_EA_output_file = output_file[:(string.rfind(output_file, "."))] + "_CM_EA.csv"

#print sys.argv

###############################################################################
#Read in Cartesian coordinates
###############################################################################

#Open input file
co_file = file(input_file, "r")

#Atom Position Vector Arrays
H1 = []
H2 = []
Oxy = []

#Ignore first two lines of the file
junk = co_file.readline()
junk = co_file.readline()

testline = co_file.readline()

while testline != "":
    # Eliminate white space and split line
    in_str = testline.strip()
    in_str = in_str.split()

    #Convert coordinates based on atom type

    if (in_str[0] == "H" and len(H1)==len(H2)):

        #Convert coordinates to vectors for the H1's
        H1.append(Geo.Vector(float(in_str[1]),float(in_str[2]),float(in_str[3])))
        print "H1"

    elif (in_str[0] == "H" and len(H2)==(len(H1)-1)):

        #Convert coordinates to vectors for the H2's
        H2.append(Geo.Vector(float(in_str[1]),float(in_str[2]),float(in_str[3])))
        print "H2"

    elif (in_str[0] == "O"):

        #Convert coordinates to vectors for the O's
        Oxy.append(Geo.Vector(float(in_str[1]),float(in_str[2]),float(in_str[3])))
        print "O"

    else:

        raise CorrelationError("Can't parse file.")
        

    #Read the next line for the while loop test condition
    testline = co_file.readline()

co_file.close()

print "Atom positions read from file"
#print Oxy

###############################################################################
#Determine the atom correlations for the water molecules
#For clathrates, all the H's are within 1A of the oxygen atom within 
#	the same water molecule, but further away from any other 
#	oxygen atom. Thus, this can be used to separate and sort out
#	which atoms belong in which molecule. The index number of the 
#	oxygen atom is taken as the index number of the molecule.
###############################################################################

#Use 2D lists to store the vectors
OH1_all = [[]]
OH2_all = [[]]

#These lists store the H1 (or H2) indices that correlate with each O
OH1_all_cor = []
OH2_all_cor = []

for waterNo in range(0,len(H1)):
    for HydNo in range(0,len(H1)):
        #Calculate vectors from O to H atoms
        OH1_all[waterNo].append(H1[HydNo]-Oxy[waterNo])
        OH2_all[waterNo].append(H2[HydNo]-Oxy[waterNo])
        #print OH1_all


        #Test if the OH1 distance is less than 1A, if so, store the index
        if (OH1_all[waterNo][HydNo].length() < 1) and (len(OH1_all_cor) < (waterNo+1)):
            OH1_all_cor.append(HydNo)

        #If there is already a value, there are too many H's correlated to the oxygen
        elif (OH1_all[waterNo][HydNo].length() < 1) and (len(OH1_all_cor) >= (waterNo+1)):
            print "Too many H1's for O%d" % (waterNo+1)
            raise CorrelationError("Too many H1's for O%d" % (waterNo+1))

        #Test if the OH2 distance is less than 1A, if so, store the index
        if (OH2_all[waterNo][HydNo].length() < 1) and (len(OH2_all_cor) < (waterNo+1)):
            OH2_all_cor.append(HydNo)
            

        #If there is already a value, there are too many H's correlated to the oxygen
        elif (OH2_all[waterNo][HydNo].length() < 1) and (len(OH2_all_cor) >= (waterNo+1)):
            print "Too many H2's for O%d" % (waterNo+1)
            raise CorrelationError("Too many H2's for O%d" % (waterNo+1))
            
    #Add empty list to end of previous list to begin a new column
    OH1_all.append([])
    OH2_all.append([])
    
    #print OH1_all[waterNo+1]
    #print OH1_all[waterNo]
    #if waterNo >= 2: print OH1_all[waterNo-1]

print "Atoms correlated"

#for waterNo in range(0,len(H1)):
#    print "O%d correlates to H1:%d H2:%d" % (waterNo+1,(OH1_all_cor[waterNo])+1,
#                                             (OH2_all_cor[waterNo])+1)

###############################################################################
#Calculate internal water molecule vectors OH1, OH2, and HH
###############################################################################
OH1 = []
OH2 = []
HH = []
for waterNo in range(0,len(H1)):
    OH1.append(H1[OH1_all_cor[waterNo]] - Oxy[waterNo])
    OH2.append(H2[OH2_all_cor[waterNo]] - Oxy[waterNo])
    HH.append(H2[OH2_all_cor[waterNo]] - H1[OH1_all_cor[waterNo]])

#print OH1
#print OH2
#print HH

print "Internal water molecule vectors calculated"

###############################################################################
#Determine whether the water has a dangling bond or not
#	by checking whether it points radially or not
#	by dotting the OH vector with the O position vector
#	and dividing by the product of the vector magnitudes
#	If the result is near 1 or -1, the vectors are nearly
#	parallel and, since the O vector points radially,
#	the OH vector points radially. If the value is near
#	zero, the vectors are perpendicular and the H is not
#	dangling.
#	In summary:
#	x_O.x_OH/(|x_O||x_OH|) = cos theta
#	if cos theata = 0, perpendicular and not dangling
#		      = 1, parallel and dangling out
#		      = -1, anti-parallel and dangling in
###############################################################################
num_dangling = [0]*len(H1)
if (sys.argv[3] == 2):
    for waterNo in range(0,len(H1)):
        #Get cos theta
        OH1_O_cos = (Oxy[waterNo]*OH1[waterNo])/(Oxy[waterNo].length()*OH1[waterNo].length())
        OH2_O_cos = (Oxy[waterNo]*OH2[waterNo])/(Oxy[waterNo].length()*OH2[waterNo].length())

        #print OH1_O_cos
        #print OH2_O_cos


        #Check to see if the H-bond is dangling, use cos theta = 0.5 as the threshold
        thresh = 0.5
        if np.absolute(OH1_O_cos) > thresh:
            num_dangling[waterNo] += 1

        if np.absolute(OH2_O_cos) > thresh:
            num_dangling[waterNo] += 1


    #print num_dangling
    print "Number of dangling H-bonds counted"

###############################################################################
#Determine the nearest neighbours of each O (the closest three)
###############################################################################

#Use a 2D list to store the vectors
OO_all = [[]]

#This 2D list stores the OO distances
OO_dist_all = [[]]

#This 2D list stores the three nearest neighbours
OO_cor = [[]]

for waterNo in range(0,len(H1)):
    for ONo in range(0,len(H1)):
        #Calculate OO vectors and distances
        OO_all[waterNo].append(Oxy[ONo]-Oxy[waterNo])
        OO_dist_all[waterNo].append(OO_all[waterNo][ONo].length())
        
            
    #Add empty list to end of previous list to begin a new column
    OO_all.append([])
    OO_dist_all.append([])

#print OO_dist_all

#Sort list and take first three interatomic distances
#This works out to be indices 1,2,3 and not 0,1,2 as the distance to itself
#is always the shortest (and is stored in index 0)
for waterNo in range(0,len(H1)):
    OO_dist_sorted = sorted(OO_dist_all[waterNo])

    #print OO_dist_sorted

    if (OO_dist_all[waterNo].count(OO_dist_sorted[1]) > 1 or
        OO_dist_all[waterNo].count(OO_dist_sorted[2]) > 1 or
        OO_dist_all[waterNo].count(OO_dist_sorted[3]) > 1):

        raise CorrelationError("Duplicate distances for O%d" % (waterNo))
    
    elif (OO_dist_all[waterNo].count(OO_dist_sorted[1]) == 0  or
        OO_dist_all[waterNo].count(OO_dist_sorted[2]) == 0 or
        OO_dist_all[waterNo].count(OO_dist_sorted[3]) == 0):

        raise CorrelationError("Distances not found for O%d" % (waterNo))

    else:
        #Store the INDICES of the nearest neighbours
        OO_cor[waterNo].append(OO_dist_all[waterNo].index(OO_dist_sorted[1]))
        OO_cor[waterNo].append(OO_dist_all[waterNo].index(OO_dist_sorted[2]))
        OO_cor[waterNo].append(OO_dist_all[waterNo].index(OO_dist_sorted[3]))

    #Add empty list to end of previous list to begin a new column
    OO_cor.append([])
    
print "Oxygen nearest neighbours determined"

###############################################################################
#Calculate the distances and angles between the nearest neighbours
###############################################################################
OO_vec = [[]]
OO_dist = [[]]

#Use a 2D-list to store the angles, as follows:
#Index:     Angle:
#  0        O1-O-O2
#  1        O2-O-O3
#  2        O3-O-O1

OOO_angle = [[]]

#Calculate the vectors, angles, and distances between each O appropriately
for waterNo in range(0,len(H1)):
    #Vectors
    OO_vec[waterNo].append(Oxy[OO_cor[waterNo][0]] - Oxy[waterNo])
    OO_vec[waterNo].append(Oxy[OO_cor[waterNo][1]] - Oxy[waterNo])
    OO_vec[waterNo].append(Oxy[OO_cor[waterNo][2]] - Oxy[waterNo])

    #Distances
    OO_dist[waterNo].append(OO_vec[waterNo][0].length())
    OO_dist[waterNo].append(OO_vec[waterNo][1].length())
    OO_dist[waterNo].append(OO_vec[waterNo][2].length())

    #Angles
    OOO_angle[waterNo].append(OO_vec[waterNo][0].angle(OO_vec[waterNo][1]))
    OOO_angle[waterNo].append(OO_vec[waterNo][1].angle(OO_vec[waterNo][2]))
    OOO_angle[waterNo].append(OO_vec[waterNo][2].angle(OO_vec[waterNo][0]))

    #Append empty list to add another column
    OO_vec.append([])
    OO_dist.append([])
    OOO_angle.append([])

print "Oxygen-Oxygen geometries calculated"

###############################################################################
#Determine toward which O the H's point and the angle from the O-O vector
###############################################################################

#Note: a value of -1 means the H is dangling and a non-negative value is the index
# of the oxygen
H1_point = []
H1_angle = []
H2_point = []
H2_angle = []

#Compare angles and then store the index of the oxygen to which the H points
for waterNo in range(0,len(H1)):
    #Get pointing direction and angle for H2
    H1O1_angle_deg = np.degrees(OH1[waterNo].angle(OO_vec[waterNo][0]))
    H1O2_angle_deg = np.degrees(OH1[waterNo].angle(OO_vec[waterNo][1]))
    H1O3_angle_deg = np.degrees(OH1[waterNo].angle(OO_vec[waterNo][2]))
    
    if ((H1O1_angle_deg <= 1)
        and (H1O2_angle_deg > 1)
        and (H1O3_angle_deg > 1)):
        H1_point.append(OO_cor[waterNo][0])
        H1_angle.append(H1O1_angle_deg)
        
    elif ((H1O1_angle_deg > 1)
        and (H1O2_angle_deg <= 1)
        and (H1O3_angle_deg > 1)):
        H1_point.append(OO_cor[waterNo][1])
        H1_angle.append(H1O2_angle_deg)
        
    elif ((H1O1_angle_deg > 1)
        and (H1O2_angle_deg > 1)
        and (H1O3_angle_deg <= 1)):
        H1_point.append(OO_cor[waterNo][2])
        H1_angle.append(H1O3_angle_deg)
        
    elif ((H1O1_angle_deg > 1)
        and (H1O2_angle_deg > 1)
        and (H1O3_angle_deg > 1)):
        H1_point.append(-1)
        H1_angle.append(0)
        
    else:
        raise CorrelationError("Can't determine direction of OH1 bond for O%d" % (waterNo))

    #Get pointing direction and angle for H2
    H2O1_angle_deg = np.degrees(OH2[waterNo].angle(OO_vec[waterNo][0]))
    H2O2_angle_deg = np.degrees(OH2[waterNo].angle(OO_vec[waterNo][1]))
    H2O3_angle_deg = np.degrees(OH2[waterNo].angle(OO_vec[waterNo][2]))
    
    if ((H2O1_angle_deg <= 1)
        and (H2O2_angle_deg > 1)
        and (H2O3_angle_deg > 1)):
        H2_point.append(OO_cor[waterNo][0])
        H2_angle.append(H2O1_angle_deg)
        
    elif ((H2O1_angle_deg > 1)
        and (H2O2_angle_deg <= 1)
        and (H2O3_angle_deg > 1)):
        H2_point.append(OO_cor[waterNo][1])
        H2_angle.append(H2O2_angle_deg)
        
    elif ((H2O1_angle_deg > 1)
        and (H2O2_angle_deg > 1)
        and (H2O3_angle_deg <= 1)):
        H2_point.append(OO_cor[waterNo][2])
        H2_angle.append(H2O3_angle_deg)
        
    elif ((H2O1_angle_deg > 1)
        and (H2O2_angle_deg > 1)
        and (H2O3_angle_deg > 1)):
        H2_point.append(-1)
        H2_angle.append(0)
        
    else:
        raise CorrelationError("Can't determine direction of OH2 bond for O%d" % (waterNo))

print "OH bond directions and angles determined"

###############################################################################
#Calculate the centre of mass for each water molecule
###############################################################################
m_tot = m_HYD1 + m_HYD2 + m_OXY
m_centre = []

for waterNo in range(0,len(H1)):
    m_vec = H1[OH1_all_cor[waterNo]]*m_HYD1 + H2[OH2_all_cor[waterNo]]*m_HYD2 + Oxy[waterNo]*m_OXY
    m_vec = m_vec/m_tot
    m_centre.append(m_vec)


###############################################################################
#Calculate the Euler angles for each water molecule
###############################################################################
# Take the water molecule's reference frame as follows:
# z-axis points from the centre of mass (CM) to the O atom
# x-axis points such that CM-H2*x-axis > 0
# y-axis = z-axis X x-axis (thus also CM-H2 X CM-H1)
#Thus,   O       Z
#       / \      ^
#      H1 H2     |->X
# i.e., the water molecule sits on the ZX plane with the CM as the origin
#       and O is in the positive Z direction
#       and H2 is in quadrant 4 of the ZX plane
#So, calculate Z as Z_ref = (CM-O)/|CM-O|
#              Y as Y_ref = (CM-H2 X CM-H1)/|CM-H2 X CM-H1|
#              X as X_ref = Y_ref X Z_ref
#Then, if the Euler angles are defined as:
# alpha = precession
# beta  = nutation
# gamma = intrinsic rotation
# They can be found as follows:
# (taken from http://en.wikipedia.org/wiki/Euler_angles#Angles_of_a_given_frame:_Geometric_derivation)
# cos(alpha) = -Z_2 / sqrt(1-(Z_3)^2)
# cos(beta)  = Z_3
# cos(gamma) = Y_3 / sqrt(1-(Z_3)^2),
#     where Z = (Z_1, Z_2, Z_3) and Y = (Y_1, Y_2, Y_3)
# Note: numpy.arccos(x) is defined as returning a value in [0, pi], however both
#       alpha and gamma are defined from [0, 2pi]. Thus, use numpy.arctan2(x1,x2).
#       From the same Wikipedia reference, alpha and gamma are thus defined as
#       follows:
# alpha =  arctan2(Z_1,-Z_2)
# gamma = arctan2(X_3,Y_3)
# Then, in MMTK, the molecule position and orientation can be entered with:
# postion as centre of mass vector
# orientation as the the composition of the extrinsic rotations: zxz,
#   using (gamma, beta, alpha) as the corresponding angles
CM_O = []
CM_H1 = []
CM_H2 = []

alpha_W = []
beta_W = []
gamma_W = []

for waterNo in range(0, len(H1)):
    #Calculate vectors from CM to each atom in the water molecule
    CM_O.append(Oxy[waterNo]-m_centre[waterNo])
    CM_H1.append(H1[OH1_all_cor[waterNo]]-m_centre[waterNo])
    CM_H2.append(H2[OH2_all_cor[waterNo]]-m_centre[waterNo])

    #Calculate Z and Y for the molecule reference frame
    Z_ref = CM_O[waterNo] / CM_O[waterNo].length()
    temp = CM_H2[waterNo].cross(CM_H1[waterNo])
    Y_ref = temp/temp.length()
    X_ref = Y_ref.cross(Z_ref)

    #Calculate Euler angles
    beta_W.append(np.arccos(Z_ref[2]))
    
    alpha = np.arctan2(Z_ref[0], -Z_ref[1])
    #convert to range [0, 2*pi)
    if (alpha < 0):
        alpha_W.append(alpha+2*np.pi)
    else:
        alpha_W.append(alpha)

        
    gamma = np.arctan2(X_ref[2], Y_ref[2])
    #convert to range [0, 2*pi)
    if (gamma < 0):
        gamma_W.append(gamma+2*np.pi)
    else:
        gamma_W.append(gamma)
    
    
    
###############################################################################
#Print the geometrical data to a .csv file for manipulation in Excel
###############################################################################

#Open the .csv file and replace any existing files
data_file = open(output_file,'w')

#####
#Write out the data for the water molecule geometries
#####

#Write the Header
data_file.write("Bond Lengths and Angles for the Clathrate from file " + input_file + "\n")
data_file.write("Molecule #, O-H1 (A), O-H2 (A), H1-H2 (A), H-O-H Angle (deg),"
                + " Number of Dangling Bonds, H1 #, H2 #" + "\n")

#Write out the data
for waterNo in range(0,len(H1)):
    #Molecule Number
    data_file.write("%d, " % (waterNo+1))

    #OH1 Vector length
    data_file.write("%.15G, " % (OH1[waterNo].length()))

    #OH2 Vector length
    data_file.write("%.15G, " % (OH2[waterNo].length()))

    #HH Vector length
    data_file.write("%.15G, " % (HH[waterNo].length()))

    #HOH Angle
    HOHangle_rad = OH2[waterNo].angle(OH1[waterNo])
    HOHangle_deg = np.degrees(HOHangle_rad)
    data_file.write("%.15G, " % (HOHangle_deg))

    #Number of Dangling H-bonds
    data_file.write("%d, " % (num_dangling[waterNo]))
    
    #H1 Number
    data_file.write("%d, " % ((OH1_all_cor[waterNo])+1))

    #H2 Number
    data_file.write("%d, " % ((OH2_all_cor[waterNo])+1))

    #New Line Character
    data_file.write("\n")

#####
#Write out the data for the O-O geometries
#####

#Write the Header
data_file.write("\n")
data_file.write("Oxygen separation and Angles for the Clathrate from file " + input_file + "\n")
data_file.write("Molecule #, O1 #, O2 #, O3 #, O-O1 (A), O-O2 (A), O-O3 (A), "
                + "O1-O-O2 Angle (deg), O2-O-O3 Angle (deg), O3-O-O1 Angle (deg), "
                + "OH1 Direction, OH1:OO Angle (deg), OH2 Direction, OH2:OO Angle (deg)" + "\n")

#Write out the data
for waterNo in range(0,len(H1)):
    #Molecule Number
    data_file.write("%d, " % (waterNo+1))

    #O1 Number
    data_file.write("%d, " % (OO_cor[waterNo][0]+1))

    #O2 Number
    data_file.write("%d, " % (OO_cor[waterNo][1]+1))

    #O3 Number
    data_file.write("%d, " % (OO_cor[waterNo][2]+1))

    #OO1 Vector length
    data_file.write("%.15G, " % (OO_dist[waterNo][0]))

    #OO2 Vector length
    data_file.write("%.15G, " % (OO_dist[waterNo][1]))

    #OO3 Vector length
    data_file.write("%.15G, " % (OO_dist[waterNo][2]))

    #O1-O-O2 Angle
    data_file.write("%.15G, " % (np.degrees(OOO_angle[waterNo][0])))

    #O2-O-O3 Angle
    data_file.write("%.15G, " % (np.degrees(OOO_angle[waterNo][1])))

    #O3-O-O1 Angle
    data_file.write("%.15G, " % (np.degrees(OOO_angle[waterNo][2])))

    #OH1 Direction (to which O the vector points)
    data_file.write("%d, " % (H1_point[waterNo]+1))

    #OH1:OO Angle
    data_file.write("%.15G, " % (H1_angle[waterNo]))

    #OH2 Direction (to which O the vector points)
    data_file.write("%d, " % (H2_point[waterNo]+1))

    #OH2:OO Angle
    data_file.write("%.15G, " % (H2_angle[waterNo]))
    
    #New Line Character
    data_file.write("\n")

data_file.write("Note: the OH# Direction is the oxygen atom to which the OH# vector points" + "\n"
                + " and a value of 0 means the bond is dangling." + "\n"
                + " i.e. if OH1 Direction is 3 then the geometry is O-H1-O3" + "\n")
data_file.write("Note: the OH#:OO angle is the angle between the OH# vector and the OO vector." + "\n"
                + " i.e. if OH1:OO is 2 deg and OH1 Direction is 3 the angle H1-O-O3" + "\n"
                + " is 2 deg. This is a measure of the H-bond linearity; the closer to 0 deg" + "\n"
                + " the more linear." + "\n")

#Close the file
data_file.close()

print ("Geometrical data written to " + output_file)

###############################################################################
#Print the centre of mass data to a .xyz file for input into VMD
###############################################################################

#Open the .xyz file and replace any existing files
CM_data_file = open(CM_output_file,'w')

#Write the header, which is the number of entries and who generated the file
CM_data_file.write("%d" % (len(H1)) + "\n")
CM_data_file.write("File generated by: " + os.path.basename(__file__) + "\n")

#Print the waterNo as W# and then the x, y, z coordinates of each water molecule
for waterNo in range(0, len(H1)):
    x_cor = m_centre[waterNo][0]
    y_cor = m_centre[waterNo][1]
    z_cor = m_centre[waterNo][2]
    CM_data_file.write("W%d %.15G %.15G %.15G" % (waterNo+1, x_cor, y_cor, z_cor) + "\n")

CM_data_file.close()

print ("Centre of Mass data written to " + CM_output_file)

###############################################################################
#Print the centre of mass and Euler angle data to a .csv file
###############################################################################

#Open the .csv file and replace any existing files
CM_EA_data_file = open(CM_EA_output_file,'w')

#Write the header, which is the number of entries and who generated the file
CM_EA_data_file.write("Centre of Mass positions and Euler Angles for file: "
                      + CM_EA_output_file + "\n")
CM_EA_data_file.write("Molecule #, CM-x (A), CM-y (A), CM-z (A), "
                      + "Precession: alpha (rad), Nutation: beta (rad), "
                      + "Instrinsic Rotation: gamma (rad)" + "\n")

#Print molecule # and then the CM coordinates and Euler Angles
#Note:
# alpha = precession
# beta  = nutation
# gamma = intrinsic rotation
for waterNo in range(0, len(H1)):
    x_cor = m_centre[waterNo][0]
    y_cor = m_centre[waterNo][1]
    z_cor = m_centre[waterNo][2]
    alpha = alpha_W[waterNo]
    beta = beta_W[waterNo]
    gamma = gamma_W[waterNo]
    CM_EA_data_file.write("%d, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G"
                       % (waterNo+1, x_cor, y_cor, z_cor, alpha, beta, gamma) + "\n")

CM_EA_data_file.close()

print ("Centre of Mass and Euler Angle data written to " + CM_EA_output_file)


