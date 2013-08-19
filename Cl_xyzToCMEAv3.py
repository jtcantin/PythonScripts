###############################################################################
#Program to calculate the Centres of Mass and Euler Angles of a Clathrate Hydrate
###############################################################################
#Program created on 12 Apr 2013 by Joshua Cantin and based on "ClGeoCalcv2.py"
# Version 1 completed on 16 Apr 2013 by Joshua Cantin
# Version 2 completed on 17 Apr 2013 by Joshua Cantin
#           - Shifted reading in the Cartesian Coordinates and Atom Correlation to a module
# Version 3 completed on 17 May 2013 by Joshua Cantin
#           - Added an xyz file for the molecules, but with no beads
###############################################################################

###############################################################################
# To run: pydev Cl_xyzToCMEAv2.py [input_file] [output_file] [N_trans] [N_rot]
# [input_file]:  An xyz file containing the hydrogens and oxygens of each water molecule.
# [output_file]: An xyz file to contain the centre of mass coordinates.
# [N_trans]:     The number of translational beads; must be an integer multiple of N_rot.
# [N_rot]:       The number of rotational beads; must be a factor of N_trans.
#################
#Output: Five files
# 1. An xyz file containing the centre of mass (CM) coordinates of each water molecule.
#    Ends with _CM.xyz
# 2. A csv file containing the CM and Euler Angles of each water molecule (zxz convention).
#    Ends with _CM_EA_zxz.csv
# 3. A csv file containing the CM and Euler Angles of each water molecule (ZY'Z'' convention).
#    Ends with _CM_EA_ZYZ.csv
# 4. An xyz file containing the CM and Euler Angles of each bead of every water molecule (ZY'Z'' convention).
#    Ends with _BEADS_CM_EA_ZYZ.xyz
# 5. An xyz file containing the CM and Euler Angles of every water molecule (ZY'Z'' convention).
#    Ends with _CM_EA_ZYZ.xyz
###############################################################################

###############################################################################
#Assumptions that this code makes:
#1. That the OH bond is less than 1A and that no H's are closer to the O
#	then those that are part of the water molecule
#	-> This affects the determination of which H's are part of the same
#		water molecule as the O-atom; if violated, the distances and
#		other properties calculated will be invalid
###############################################################################

#NOTE: Floating point numbers are usually represented as double in C, unless
#      the machine cannot handle double precision
#      See section 5.4 of http://docs.python.org/2/library/stdtypes.html 

import numpy as np
import Scientific as Sci
import Scientific.Geometry as Geo
import sys
import string
import os
import Clath_xyz
from MMTK import *

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

## N_a = 6.022141E23 #/mol Taken from:CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL
##                   #     PHYSICAL CONSTANTS: 2010. In CRC Handbook of Chemistry
##                   #     and Physics, 93rd ed.; CRC Press, 2012-2013.

N_a = Units.Nav

m_HYD1 = M_HYD1/N_a/1000 #kg
m_HYD2 = M_HYD2/N_a/1000 #kg
m_OXY = M_OXY/N_a/1000 #kg

print "Program Executing"

input_file = sys.argv[1] 
output_file = sys.argv[2]

N_trans = int(sys.argv[3])
N_rot = int(sys.argv[4])

if (N_trans/N_rot % 1 != 0):
    raise CorrelationError("The number of translational beads must be a multiple of the number of rotaitonal beads")

#Strip output file extension and then replace it with _CM.xyz
CM_output_file = output_file[:(string.rfind(output_file, "."))] + "_CM.xyz"
CM_EA_zxz_output_file = output_file[:(string.rfind(output_file, "."))] + "_CM_EA_zxz.csv"
CM_EA_ZYZ_output_file = output_file[:(string.rfind(output_file, "."))] + "_CM_EA_ZYZ.csv"
CM_EA_ZYZ_xyz_output_file = output_file[:(string.rfind(output_file, "."))] + "_BEADS_CM_EA_ZYZ.xyz"
CM_EA_ZYZ_xyz_normal_output_file = output_file[:(string.rfind(output_file, "."))] + "_CM_EA_ZYZ.xyz"

###############################################################################
#Read in Cartesian coordinates
###############################################################################

out_dict = Clath_xyz.readWaterAtoms(input_file)

Oxy = out_dict["O"]
Hyd = out_dict["H"]

print "Atom positions read from file."

###############################################################################
#Determine the atom correlations for the water molecules
#For clathrates, all the H's are within 1A of the oxygen atom within 
#	the same water molecule, but further away from any other 
#	oxygen atom. Thus, this can be used to separate and sort out
#	which atoms belong in which molecule. The index number of the 
#	oxygen atom is taken as the index number of the molecule.
###############################################################################

OH_length = 1 #Angstroms

out_dict = Clath_xyz.correlateWaters(Oxy, Hyd, OH_length)

#Use 2D lists to store the vectors
OH_all = out_dict["OH_all_vec"]

#These lists store the H1 (or H2) indices that correlate with each O
OH1_all_cor = out_dict["OH1_all_cor"]
OH2_all_cor = out_dict["OH2_all_cor"]
    
print "Atoms correlated."

###############################################################################
#Calculate internal water molecule vectors OH1, OH2, and HH
###############################################################################
OH1 = []
OH2 = []
HH = []
for waterNo in range(0,len(Oxy)):
    OH1.append(Hyd[OH1_all_cor[waterNo]] - Oxy[waterNo])
    OH2.append(Hyd[OH2_all_cor[waterNo]] - Oxy[waterNo])
    HH.append(Hyd[OH2_all_cor[waterNo]] - Hyd[OH1_all_cor[waterNo]])

#print OH1
#print OH2
#print HH

print "Internal water molecule vectors calculated."

###############################################################################
#Calculate the centre of mass for each water molecule
###############################################################################
m_tot = m_HYD1 + m_HYD2 + m_OXY
m_centre = []

for waterNo in range(0,len(Oxy)):
    m_vec = Hyd[OH1_all_cor[waterNo]]*m_HYD1 + Hyd[OH2_all_cor[waterNo]]*m_HYD2 + Oxy[waterNo]*m_OXY
    m_vec = m_vec/m_tot
    m_centre.append(m_vec)

###############################################################################
#Calculate the Euler angles for each water molecule
###############################################################################
    
##################################
#zxz intrinsic Convention
##################################
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

for waterNo in range(0, len(Oxy)):
    #Calculate vectors from CM to each atom in the water molecule
    CM_O.append(Oxy[waterNo]-m_centre[waterNo])
    CM_H1.append(Hyd[OH1_all_cor[waterNo]]-m_centre[waterNo])
    CM_H2.append(Hyd[OH2_all_cor[waterNo]]-m_centre[waterNo])

    #Calculate Z and Y for the molecule reference frame
    Z_ref = CM_O[waterNo] / CM_O[waterNo].length()
    temp = CM_H2[waterNo].cross(CM_H1[waterNo])
    Y_ref = temp/temp.length()
    X_temp = Y_ref.cross(Z_ref)
    X_ref = X_temp/X_temp.length()

    #Check that the vectors are of unit length
    if ((abs(Z_ref.length()-1)>np.spacing(1)) or (abs(X_ref.length()-1)>np.spacing(1)) or (abs(Y_ref.length()-1)>np.spacing(1))):
        raise CorrelationError("The zxz body-fixed frame generators are not of unit length for water ", waterNo)

    #Calculate Euler angles
    ## beta_W.append(np.arccos(Z_ref[2]))
    
    ## alpha = np.arctan2(Z_ref[0], -Z_ref[1])
    ## #convert to range [0, 2*pi)
    ## if (alpha < 0):
    ##     alpha_W.append(alpha+2*np.pi)
    ## else:
    ##     alpha_W.append(alpha)

        
    ## gamma = np.arctan2(X_ref[2], Y_ref[2])
    ## #convert to range [0, 2*pi)
    ## if (gamma < 0):
    ##     gamma_W.append(gamma+2*np.pi)
    ## else:
    ##     gamma_W.append(gamma)

        #########################
    if (abs(Z_ref[2] - 1) < np.spacing(1)): #if cos(beta)=Z_ref[2] is effectively 1, beta = 0
        beta_W.append(0.)
        #cos_beta_W.append(1.)
        
        #If beta = 0, alpha and gamma are not unique, only alpha+gamma is; use the convention that alpha = 0.
        alpha_W.append(0.)

        gamma_temp = np.arctan2(X_ref[1],Y_ref[1]) #Use arctan2 instead of coding the quadrant checking myself so as to follow the C and IEEE standards for 0 and Inf values
                                                   #Use X_ref[1] and Y_ref[1] as alpha and phi are zero, meaning that X_ref[1] = sin(gamma) and Y_ref[1] = cos(gamma)

        #Bring the range from [-pi,pi] to [0,2pi]
        if (gamma_temp < 0):
            gamma_W.append(gamma_temp+2.*np.pi)
        else:
            gamma_W.append(gamma_temp)
        
    elif (abs(Z_ref[2] + 1) < np.spacing(1)): #if cos(beta)=Z_ref[2] is effectively -1, beta = pi
        beta_W.append(np.pi)
        #cos_beta_W.append(-1.)

        #If beta = pi, alpha and gamma is not unique, only alpha+gamma is; use the convention that alpha = 0.
        alpha_W.append(0.)

        gamma_temp = np.arctan2(-X_ref[1],-Y_ref[1]) #Use arctan2 instead of coding the quadrant checking myself so as to follow the C and IEEE standards for 0 and Inf values
                                                 #Use X_ref[1] and Y_ref[1] as alpha is zero, meaning that X_ref[1] = -sin(gamma) and Y_ref[1] = -cos(gamma)
        
        #Bring the range from [-pi,pi] to [0,2pi]
        if (gamma_temp < 0):
            gamma_W.append(gamma_temp+2.*np.pi)
        else:
            gamma_W.append(gamma_temp)
        
    else: #otherwise, beta = cos^-1(Z_ref[2])    
        beta_W.append(np.arccos(Z_ref[2])) #Returns a value between [0,pi], so no further modifications needed
        #cos_beta_W.append(Z_ref[2])

        #Calculate alpha using alpha = arctan2(rotmat(1,3), -rotmat(2,3)) = acrtan2(Z_ref[0], -Z_ref[1]) 
        alpha_temp = np.arctan2(Z_ref[0],-Z_ref[1]) #Use arctan2 instead of coding the quadrant checking myself so as to follow the C and IEEE standards for 0 and Inf values
                                                 
        #Bring the range from [-pi,pi] to [0,2pi]
        if (alpha_temp < 0):
            alpha_W.append(alpha_temp+2.*np.pi)
        else:
            alpha_W.append(alpha_temp)

        #Calculate gamma using gamma = arctan2(rotmat(3,1), rotmat(3,2)) = acrtan2(X_ref[2], Y_ref[2])
        gamma_temp = np.arctan2(X_ref[2], Y_ref[2]) #Use arctan2 instead of coding the quadrant checking myself so as to follow the C and IEEE standards for 0 and Inf values
                                                 
        #Bring the range from [-pi,pi] to [0,2pi]
        if (gamma_temp < 0):
            gamma_W.append(gamma_temp+2.*np.pi)
        else:
            gamma_W.append(gamma_temp)

##################################
#ZY'Z'' Extrinsic Convention
##################################
# Take the water molecule's reference frame as follows:
# z-axis points from the centre of mass (CM) to the O atom
# x-axis points from H1 to H2
# y-axis = z-axis X x-axis (thus also CM-H2 X CM-H1)
#Thus,   O       Z
#       / \      ^
#      H1 H2     |->X
# i.e., the water molecule sits on the ZX plane with the CM as the origin
#       and O is in the positive Z direction
#       and H2 is in quadrant 4 of the ZX plane
#So, calculate Z as Z_ref = (CM-O)/|CM-O|
#              X as X_ref = (H2-H1)/|H2-H1|
#              Y as Y_ref = Z_ref x X_ref
#Then, if the Euler angles are defined as:
# phi = precession
# theta  = nutation
# chi = intrinsic rotation
# They can be found as follows:
# (derived similarly to http://en.wikipedia.org/wiki/Euler_angles#Angles_of_a_given_frame:_Geometric_derivation)
# cos(phi) = Z_1 / sqrt(1-(Z_3)^2)
# cos(theta) = Z_3
# cos(chi) = X_3 / sqrt(1-(Z_3)^2),
#     where Z = (Z_1, Z_2, Z_3) and X = (X_1, X_2, X_3)
# Note: numpy.arccos(x) is defined as returning a value in [0, pi], however both
#       phi and chi are defined from [0, 2pi]. Thus, use numpy.arctan2(x1,x2).
#       Similarly, alpha and gamma are thus defined as
#       follows:
# phi =  arctan2(Z_2,Z_1)
# chi = arctan2(Y_3,-X_3) (equivalent to rotmat(3,2)/rotmat(3,1), where rotmat takes the body fixed frame to the space fixed frame)
# Then, in MMTK, the molecule position and orientation can be entered with:
# postion as centre of mass vector
# orientation as the the composition of the intrinsic rotations: ZY'Z'',
#   using (phi, theta, chi) as the corresponding angles
        
CM_O = []
CM_H1 = []
CM_H2 = []
H1_H2 = []

phi_W = []
theta_W = []
chi_W = []
cos_theta_W = []

for waterNo in range(0, len(Oxy)):
    #Calculate vectors from CM to each atom in the water molecule
    CM_O.append(Oxy[waterNo]-m_centre[waterNo])
    CM_H1.append(Hyd[OH1_all_cor[waterNo]]-m_centre[waterNo])
    CM_H2.append(Hyd[OH2_all_cor[waterNo]]-m_centre[waterNo])
    H1_H2.append(Hyd[OH2_all_cor[waterNo]]-Hyd[OH1_all_cor[waterNo]])

    ## #Calculate Z and Y for the molecule reference frame
    ## Z_ref = CM_O[waterNo] / CM_O[waterNo].length()
    ## X_ref = H1_H2[waterNo] / H1_H2[waterNo].length()
    ## Y_ref_temp = Z_ref.cross(X_ref) #Calculate in two steps to reduce errors in cross-product calculation
    ## Y_ref = Y_ref_temp/Y_ref_temp.length()

    # Use the below method and not the above as this gives a smaller RMS difference. Toby mentioned that this is because of numerical errors
          # modifying the OH bond lengths slightly and causing H2_H1 to not be actually parallel to the desired x-axis, while the below method
          # seems to preserve the parallelism.
    #Calculate Z and Y for the molecule reference frame
    Z_ref = CM_O[waterNo] / CM_O[waterNo].length()
    temp = CM_H2[waterNo].cross(CM_H1[waterNo])
    Y_ref = temp/temp.length() #Calculate in two steps to reduce errors in cross-product calculation
    X_temp = Y_ref.cross(Z_ref)
    X_ref = X_temp/X_temp.length() #Calculate in two steps to reduce errors in cross-product calculation

    #Check that the vectors are of unit length
    if ((abs(Z_ref.length()-1)>np.spacing(1)) or (abs(X_ref.length()-1)>np.spacing(1)) or (abs(Y_ref.length()-1)>np.spacing(1))):
        raise CorrelationError("The ZY'Z'' body-fixed frame generators are not of unit length for water ", waterNo)

    #Calculate Euler angles

    if (abs(Z_ref[2] - 1) < np.spacing(1)): #if cos(theta)=Z_ref[2] is effectively 1, theta = 0
        theta_W.append(0.)
        cos_theta_W.append(1.)
        
        #If theta = 0, phi and chi are not unique, only phi+chi is; use the convention that phi = 0.
        phi_W.append(0.)

        chi_temp = np.arctan2(X_ref[1],Y_ref[1]) #Use arctan2 instead of coding the quadrant checking myself so as to follow the C and IEEE standards for 0 and Inf values
                                                 #Use X_ref[1] and Y_ref[1] as phi is zero, meaning that X_ref[1] = sin(chi) and Y_ref[1] = cos(chi)

        #Bring the range from [-pi,pi] to [0,2pi]
        if (chi_temp < 0):
            chi_W.append(chi_temp+2.*np.pi)
        else:
            chi_W.append(chi_temp)
        
    elif (abs(Z_ref[2] + 1) < np.spacing(1)): #if cos(theta)=Z_ref[2] is effectively -1, theta = pi
        theta_W.append(np.pi)
        cos_theta_W.append(-1.)

        #If theta = pi, phi and chi is not unique, only phi+chi is; use the convention that phi = 0.
        phi_W.append(0.)

        chi_temp = np.arctan2(X_ref[1],Y_ref[1]) #Use arctan2 instead of coding the quadrant checking myself so as to follow the C and IEEE standards for 0 and Inf values
                                                 #Use X_ref[1] and Y_ref[1] as phi is zero, meaning that X_ref[1] = sin(chi) and Y_ref[1] = cos(chi)
        
        #Bring the range from [-pi,pi] to [0,2pi]
        if (chi_temp < 0):
            chi_W.append(chi_temp+2.*np.pi)
        else:
            chi_W.append(chi_temp)
        
    else: #otherwise, theta = cos^-1(Z_ref[2])    
        theta_W.append(np.arccos(Z_ref[2])) #Returns a value between [0,pi], so no further modifications needed
        cos_theta_W.append(Z_ref[2])

        #Calculate phi using phi = arctan2(rotmat(2,3), rotmat(1,3)) = acrtan2(Z_ref[1], Z_ref[0]) 
        phi_temp = np.arctan2(Z_ref[1],Z_ref[0]) #Use arctan2 instead of coding the quadrant checking myself so as to follow the C and IEEE standards for 0 and Inf values
                                                 
        #Bring the range from [-pi,pi] to [0,2pi]
        if (phi_temp < 0):
            phi_W.append(phi_temp+2.*np.pi)
        else:
            phi_W.append(phi_temp)

        #Calculate chi using chi = arctan2(rotmat(3,2), -rotmat(3,1)) = acrtan2(Y_ref[2], -X_ref[2])
        chi_temp = np.arctan2(Y_ref[2], -X_ref[2]) #Use arctan2 instead of coding the quadrant checking myself so as to follow the C and IEEE standards for 0 and Inf values
                                                 
        #Bring the range from [-pi,pi] to [0,2pi]
        if (chi_temp < 0):
            chi_W.append(chi_temp+2.*np.pi)
        else:
            chi_W.append(chi_temp)
    
###############################################################################
#Print the centre of mass data to a .xyz file for input into VMD
###############################################################################

#Open the .xyz file and replace any existing files
CM_data_file = open(CM_output_file,'w')

#Write the header, which is the number of entries and who generated the file
CM_data_file.write("%d" % (len(Oxy)) + "\n")
CM_data_file.write("File generated by: " + os.path.basename(__file__) + "\n")

#Print the waterNo as W# and then the x, y, z coordinates of each water molecule
for waterNo in range(0, len(Oxy)):
    x_cor = m_centre[waterNo][0]
    y_cor = m_centre[waterNo][1]
    z_cor = m_centre[waterNo][2]
    CM_data_file.write("W%d %.5e %.5e %.5e" % (waterNo+1, x_cor, y_cor, z_cor) + "\n")

CM_data_file.close()

print ("Centre of Mass data written to " + CM_output_file)

###############################################################################
#Print the centre of mass and Euler angle (zxz) data to a .csv file
###############################################################################

#Open the .csv file and replace any existing files
CM_EA_data_file = open(CM_EA_zxz_output_file,'w')

#Write the header, which is the number of entries and who generated the file
CM_EA_data_file.write("Centre of Mass positions and Euler Angles with the zxz convention for file: "
                      + CM_EA_zxz_output_file + "\n")
CM_EA_data_file.write("Molecule #, CM-x (A), CM-y (A), CM-z (A), "
                      + "Precession: alpha (rad), Nutation: beta (rad), "
                      + "Instrinsic Rotation: gamma (rad)" + "\n")

#Print molecule # and then the CM coordinates and Euler Angles
#Note:
# alpha = precession
# beta  = nutation
# gamma = intrinsic rotation
for waterNo in range(0, len(Oxy)):
    x_cor = m_centre[waterNo][0]
    y_cor = m_centre[waterNo][1]
    z_cor = m_centre[waterNo][2]
    alpha = alpha_W[waterNo]
    beta = beta_W[waterNo]
    gamma = gamma_W[waterNo]
    CM_EA_data_file.write("%d, %.5e, %.5e, %.5e, %.5e, %.5e, %.5e"
                       % (waterNo+1, x_cor, y_cor, z_cor, alpha, beta, gamma) + "\n")

CM_EA_data_file.close()

print ("Centre of Mass and Euler Angle (zxz) data written to " + CM_EA_zxz_output_file)

###############################################################################
#Print the centre of mass and Euler angle (ZY'Z'') data to a .csv file
###############################################################################

#Open the .csv file and replace any existing files
CM_EA_data_file = open(CM_EA_ZYZ_output_file,'w')

#Write the header, which is the number of entries and who generated the file
CM_EA_data_file.write("Centre of Mass positions and Euler Angles with the ZY'Z'' convention for file: "
                      + CM_EA_ZYZ_output_file + "\n")
CM_EA_data_file.write("Molecule #, CM-x (A), CM-y (A), CM-z (A), "
                      + "Precession: phi (rad), Nutation: theta (rad), "
                      + "Instrinsic Rotation: chi (rad)" + "\n")

#Print molecule # and then the CM coordinates and Euler Angles
#Note:
# phi = precession
# theta  = nutation
# chi = intrinsic rotation
for waterNo in range(0, len(Oxy)):
    x_cor = m_centre[waterNo][0]
    y_cor = m_centre[waterNo][1]
    z_cor = m_centre[waterNo][2]
    phi = phi_W[waterNo]
    theta = theta_W[waterNo]
    chi = chi_W[waterNo]
    CM_EA_data_file.write("%d, %.5e, %.5e, %.5e, %.5e, %.5e, %.5e"
                       % (waterNo+1, x_cor, y_cor, z_cor, phi, theta, chi) + "\n")

CM_EA_data_file.close()

print ("Centre of Mass and Euler Angle (ZY'Z'') data written to " + CM_EA_ZYZ_output_file)

###############################################################################
#Print the centre of mass and Euler angle (ZY'Z'') data to an .xyz file with translational and rotational beads
###############################################################################

#Open the .xyz file and replace any existing files
CM_EA_data_file = open(CM_EA_ZYZ_xyz_output_file,'w')

#Write the header, which is the number of entries and who generated the file
CM_EA_data_file.write("%d" % (N_trans*len(Oxy)) + "\n")
CM_EA_data_file.write("# Format: [name_id]; CM-x (A); phi (rad); CM-y (A); cos(theta); CM-z (A); chi (rad); Euler angles for ZY'Z'' convention" + "\n")

#Print molecule # and then the CM coordinates and Euler Angles
#Note:
# phi = precession
# theta = nutation
# chi = intrinsic rotation
for waterNo in range(0, len(Oxy)):
    for transBead in range(0, N_trans):

        #CM coordinates
        x_cor = m_centre[waterNo][0]
        y_cor = m_centre[waterNo][1]
        z_cor = m_centre[waterNo][2]
        
        if (transBead < N_rot):
            phi = phi_W[waterNo]
            cos_theta = cos_theta_W[waterNo]
            chi = chi_W[waterNo]
        else:
            phi = 0.
            cos_theta = 1.
            chi = 0.

        CM_EA_data_file.write("H2O%d %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e"
                       % (waterNo+1, x_cor, phi, y_cor, cos_theta, z_cor, chi) + "\n")

CM_EA_data_file.close()

print ("Centre of Mass and Euler Angle (ZY'Z'') data written to " + CM_EA_ZYZ_xyz_output_file)


###############################################################################
#Print the centre of mass and Euler angle (ZY'Z'') data to an .xyz file
###############################################################################

#Open the .xyz file and replace any existing files
CM_EA_data_file = open(CM_EA_ZYZ_xyz_normal_output_file,'w')

#Write the header, which is the number of entries and who generated the file
CM_EA_data_file.write("%d" % (len(Oxy)) + "\n")
CM_EA_data_file.write("# Format: [name_id]; CM-x (A); phi (rad); CM-y (A); theta (rad); CM-z (A); chi (rad); Euler angles for ZY'Z'' convention" + "\n")

#Print molecule # and then the CM coordinates and Euler Angles
#Note:
# phi = precession
# theta = nutation
# chi = intrinsic rotation
for waterNo in range(0, len(Oxy)):

        #CM coordinates
        x_cor = m_centre[waterNo][0]
        y_cor = m_centre[waterNo][1]
        z_cor = m_centre[waterNo][2]
        

        phi = phi_W[waterNo]
        theta = theta_W[waterNo]
        chi = chi_W[waterNo]

        CM_EA_data_file.write("H2O%d %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e"
                       % (waterNo+1, x_cor, phi, y_cor, theta, z_cor, chi) + "\n")

CM_EA_data_file.close()

print ("Centre of Mass and Euler Angle (ZY'Z'') data written to " + CM_EA_ZYZ_xyz_normal_output_file)


