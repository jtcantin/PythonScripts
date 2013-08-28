####################################################
####################################################
# Usage: pydev genLanczosHvParamFiles.py [outputDir]
####################################################
####################################################
# Written by Joshua Cantin 27 Aug 2013

from numpy import *
import sys
import datetime
import os
import string

####################################################
#Fixed Hv Parameters
####################################################
quadratureConvergenceStudy = "FALSE"
#geometryFilename = "sII_smallCage11_CM_EA_ZYZ.xyz"
totalMass = 2.01588 #g/mol
rotationalConstant = 0.7096487130193394 #kJ/mol
simulationTitlePrefix = "Alavi_SPCE_smallCage11"

rotationalSymmetry = "Even"

####################################################
#Miscellaneous Parameters
####################################################

nthreads = 8
executable = "LanczosTOBY"
H2Allowance = 0.08 #nm


####################################################
#Fixed Lanczos Parameters
####################################################
EigenvalueOpenLowerLimit = -10 #kJ/mol
EigenvalueClosedUpperLimit = 10 #kJ/mol
NumberOfEigenvectorsWord = "none"
NumberOfEigenvectorsNum = 3
HvCalculator = "linRotCartSph_NoQuad_Alavi_SPCE"
OutputDirectoryPath = "/scratch/jtcantin/"

####################################################
#Base Hv Parameters
####################################################

baseHvParamDict = { "gridPoints":					21,
                    "gridSize":						0.4, #nm
					"gridSpacing":					0.02, #nm
                    "lMax":							5,
                    "thetaPoints":					70,
                    "phiPoints":					70,
                    "potentialGridPoints":			241,
                    "potentialGridSize":			0.48, #nm
					"potentialGridSpacing":			0.002, #nm
                    "totalMass":					totalMass,
                    "rotationalConstant":			rotationalConstant,
                    "geometryFilename":				"sII_smallCage11_CM_EA_ZYZ.xyz",
                    "simulationTitle":				simulationTitlePrefix + "_baseParam",
                    "ceilingPotential":				18., #kJ/mol
                    "quadratureConvergenceStudy":	quadratureConvergenceStudy,
                    "rotationalSymmetry":			rotationalSymmetry
                    }
					

####################################################
#Base Lanczos Parameters
####################################################

baseLanczosParamDict = { "lanczosIterations":		10000,
                    "EigenvalueOpenLowerLimit":		EigenvalueOpenLowerLimit,
                    "EigenvalueClosedUpperLimit":	EigenvalueClosedUpperLimit,
                    "NumberOfEigenvectorsWord":		NumberOfEigenvectorsWord,
                    "NumberOfEigenvectorsNum":		NumberOfEigenvectorsNum,
                    "HvCalculator":					HvCalculator,
                    "OutputDirectoryPath":			OutputDirectoryPath
                    }

####################################################
#Parameter Ranges
####################################################
#Syntax: numpy.r_[first:last+step:step]
gridSizeArray = r_[0.4:0.5:0.1] #nm
gridSpacingArray = r_[0.02:0.03:0.01] #nm
lMaxArray = r_[5:6:1]
thetaPointsArray = r_[70:80:10]
phiPointsArray = r_[70:80:10]
ceilingPotentialArray = r_[18.:19.:1.] #kJ/mol
#potentialGridSizeArray = gridSizeArray + 0.08 #nm
potentialGridSpacingArray = array([0.0005, 0.001, 0.002, 0.003, 0.004])#r_[0.001:0.004:0.001] #nm
geometryFilenameArray = array(["sII_smallCage11_CM_EA_ZYZ.xyz"])
lanczosIterationsArray = r_[10000:20000:10000]

#Round arrays appropriately
gridSizeArray = around(gridSizeArray, 15)
gridSpacingArray = around(gridSpacingArray, 15)
lMaxArray = rint(lMaxArray).astype(int)
thetaPointsArray = rint(thetaPointsArray).astype(int)
phiPointsArray = rint(phiPointsArray).astype(int)
ceilingPotentialArray = around(ceilingPotentialArray, 15)
#potentialGridSizeArray = around(potentialGridSizeArray, 15)
potentialGridSpacingArray = around(potentialGridSpacingArray, 15)
lanczosIterationsArray = rint(lanczosIterationsArray).astype(int)

rangesList = [	gridSizeArray,
				gridSpacingArray,
				lMaxArray,
				thetaPointsArray,
				phiPointsArray,
				ceilingPotentialArray,
				#potentialGridSizeArray,
				potentialGridSpacingArray,
				geometryFilenameArray,
				lanczosIterationsArray
				]
				
rangesListText = [	"gridSize",
				"gridSpacing",
				"lMax",
				"thetaPoints",
				"phiPoints",
				"ceilingPotential",
				#"potentialGridSize",
				"potentialGridSpacing",
				"geometryFilename",
				"lanczosIterations"
				]

print ""
print "The parameter values are: "
print ""

print "gridSizeArray: "
print "     ", gridSizeArray
print ""

print "gridSpacingArray: "
print "     ", gridSpacingArray
print ""

#print "gridPointsArray: "
#print gridPointsArray
print "lMaxArray: "
print "     ", lMaxArray
print ""

print "thetaPointsArray: "
print "     ", thetaPointsArray
print ""

print "phiPointsArray: "
print "     ", phiPointsArray
print ""

print "ceilingPotentialArray: "
print "     ", ceilingPotentialArray
print ""

#print "potentialGridSizeArray: "
#print potentialGridSizeArray
print "potentialGridSpacingArray: "
print "     ", potentialGridSpacingArray
print ""

#print "potentialGridPointsArray: "
#print potentialGridPointsArray
print "geometryFilenameArray: "
print "     ", geometryFilenameArray
print ""

print "lanczosIterationsArray: "
print "     ", lanczosIterationsArray
print ""

print "quadratureConvergenceStudy: "
print "     ", quadratureConvergenceStudy
print ""

print "totalMass: "
print "     ", totalMass
print ""

print "rotationalConstant: "
print "     ", rotationalConstant
print ""

print "simulationTitlePrefix: "
print "     ", simulationTitlePrefix
print ""


####################################################
#Function Definitions
####################################################

# Hv Parameter file
####################################################
def writeHvParams(params, filename):
    datafile = open(filename, 'w')

    #Write header
    datafile.write("#" + "\n")
    datafile.write("#Make sure that there is no space before the equals sign and only one after; you will get an error or weird results otherwise" + "\n")

    ###########
    #Write data
    ###########

    #Cartesian system grid
    datafile.write("nx= %d" % (params["gridPoints"]) + "\n")
    datafile.write("x_max(nm)= %.15f" % (params["gridSize"]) + "\n")

    datafile.write("ny= %d" % (params["gridPoints"]) + "\n")
    datafile.write("y_max(nm)= %.15f" % (params["gridSize"]) + "\n")

    datafile.write("nz= %d" % (params["gridPoints"]) + "\n")
    datafile.write("z_max(nm)= %.15f" % (params["gridSize"]) + "\n")

    #l_max
    datafile.write("l_max= %d" % (params["lMax"]) + "\n")

    #Angular grid
    datafile.write("thetaPoints= %d" % (params["thetaPoints"]) + "\n")
    datafile.write("phiPoints= %d" % (params["phiPoints"]) + "\n")

    #Cartesian potential grid
    datafile.write("pointPotential_nx= %d" % (params["potentialGridPoints"]) + "\n")
    datafile.write("pointPotential_x_max(nm)= %.15f" % (params["potentialGridSize"]) + "\n")

    datafile.write("pointPotential_ny= %d" % (params["potentialGridPoints"]) + "\n")
    datafile.write("pointPotential_y_max(nm)= %.15f" % (params["potentialGridSize"]) + "\n")

    datafile.write("pointPotential_nz= %d" % (params["potentialGridPoints"]) + "\n")
    datafile.write("pointPotential_z_max(nm)= %.15f" % (params["potentialGridSize"]) + "\n")

    #Linear Rotor Parameters
    datafile.write("totalMass(g/mol)= %.15f" % (params["totalMass"]) + "\n")
    datafile.write("rotationalConstant(kJ/mol)= %.15f" % (params["rotationalConstant"]) + "\n")

    #Miscellaneous Parameters
    datafile.write("geometryFilename= %s" % (params["geometryFilename"]) + "\n")

    datafile.write("simulationTitle= %s" % (params["simulationTitle"]) + "\n")

    datafile.write("ceilingPotential(kJ/mol)= %.15f" % (params["ceilingPotential"]) + "\n")

    datafile.write("quadratureConvergenceStudy(TRUE/FALSE)= %s" % (params["quadratureConvergenceStudy"]) + "\n")

    datafile.write("rotationalSymmetry(All/Even/Odd)= %s" % (params["rotationalSymmetry"]) + "\n")

    #Finished
    datafile.close()
	
# Lanczos Parameter file
####################################################
def writeLanczosParams(params, filename):
	datafile = open(filename, 'w')
	
    #Write header
	datafile.write("#This is the input file for the Boundstate Lanczos Eigenvalue and Eigenvector Calculator" + "\n")
	datafile.write("#Make sure that there is no space before the equals sign and only one after; you will get an error or weird results otherwise" + "\n")
	
	###########
    #Write data
    ###########
	
	datafile.write("numberOfIterations= %d" % (params["lanczosIterations"]) + "\n")
	datafile.write("EigenvalueOpenLowerLimit(kJ/mol)= %.15f" % (params["EigenvalueOpenLowerLimit"]) + "\n")
	datafile.write("EigenvalueClosedUpperLimit(kJ/mol)= %.15f" % (params["EigenvalueClosedUpperLimit"]) + "\n")
	datafile.write("NumberOfEigenvectors(all/partial/none)= %s" % (params["NumberOfEigenvectorsWord"]) + "\n")
	datafile.write("NumberOfEigenvectors= %d" % (params["NumberOfEigenvectorsNum"]) + "\n")
	datafile.write("HvCalculator= %s" % (params["HvCalculator"]) + "\n")
	datafile.write("OutputDirectoryPath= %s" % (params["OutputDirectoryPath"]) + "\n")
	
	#Finished
	datafile.close()
	
# Build the simulation title
####################################################
def genSimTitle(pref, HvParams, LanczosParams):

	if HvParams["rotationalSymmetry"] == "Odd":
		rotSymTitle = "ortho"
	elif HvParams["rotationalSymmetry"] == "Even":
		rotSymTitle = "para"
	else:
		print "ERROR: rotational symmetry ill-defined"
		sys.exit()

	simTitle = pref + "_%s_x%s_dx%s_l%s_t%s_p%s_px%s_pdx%s_cp%s_L%s" % (rotSymTitle, 
									string.replace(str(HvParams["gridSize"]), ".", ""),
									string.replace(str(HvParams["gridSpacing"]), ".", ""),
									str(HvParams["lMax"]),
									str(HvParams["thetaPoints"]),
									str(HvParams["phiPoints"]),
									string.replace(str(HvParams["potentialGridSize"]), ".", ""),
									string.replace(str(HvParams["potentialGridSpacing"]), ".", ""),
									string.replace(str(HvParams["ceilingPotential"]), ".", ""),
									str(LanczosParams["lanczosIterations"])
									)
	return simTitle
	
# Append to shell script
####################################################
def appendShellFile(executable, nthreads, HvParams, HvFilename, LanczosFilename, logpath, filename):
	datafile = open(filename, 'a')
	
	datafile.write("\n")
	datafile.write("export OMP_NUM_THREADS=%d" % (nthreads) + "\n")
	datafile.write("echo %s" % (HvParams["simulationTitle"]) + "\n")
	datafile.write("time ./%s %s %s" % (executable, LanczosFilename, HvFilename))
	datafile.write(" > %s/%s" % (logpath, HvParams["simulationTitle"] + "_log") + "\n")

	#Finished
	datafile.close()

####################################################
####################################################
#Main Program
####################################################
####################################################

####################################################
#Get and make output directory
####################################################
outDir = sys.argv[1] + "_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
os.system("mkdir ./" + outDir)

HvFilePrefix = "HvInputFile"
LanzcosFilePrefix = "LanczosInputFile"

#Make a directory to contain the files to be placed in Compilation
outDirComp =  outDir + "/" + "Compilation"
os.system("mkdir -p ./" + outDirComp)

#Write out base parameters files
baseHvParamDict["simulationTitle"] = genSimTitle(simulationTitlePrefix, baseHvParamDict, baseLanczosParamDict) + "_Base"

HvFileTitle = HvFilePrefix + "_" + baseHvParamDict["simulationTitle"] + ".txt"
LanczosFileTitle = LanzcosFilePrefix + "_" + baseHvParamDict["simulationTitle"] + ".txt"

writeHvParams(baseHvParamDict, outDir + "/" + HvFileTitle)
writeHvParams(baseHvParamDict, outDirComp + "/" + HvFileTitle)

writeLanczosParams(baseLanczosParamDict, outDir + "/" + LanczosFileTitle)
writeLanczosParams(baseLanczosParamDict, outDirComp + "/" + LanczosFileTitle)

appendShellFile(executable, nthreads, baseHvParamDict, HvFileTitle, LanczosFileTitle, OutputDirectoryPath, outDirComp + "/" + "Run.sh")

#Generate parameter files for each range of values, using the base parameters as the starting point for the variations
cnt = 1 #start at 1 to include base parameters

print ""
print "Writing parameter files for:"

for i in range(0,len(rangesList)) :
	
	print ""
	print rangesListText[i] + ":"
	
	#Don't make a separate calculation for a single parameter
	if len(rangesList[i]) <= 1:
		print "     Single parameter for %s, value: %s, skipping..." % (rangesListText[i], str(rangesList[i][0]))
		continue
		
	#Make a directory for the range of values
	outDirRange = outDir + "/" + rangesListText[i]
	os.system("mkdir -p ./" + outDirRange)
	
	for value in rangesList[i]:
		
		#Reset to base parameters
		HvParamsLocal = baseHvParamDict.copy()
		LanczosParamsLocal = baseLanczosParamDict.copy()
				
		if	i < (len(rangesList)-1):
			HvParamsLocal[rangesListText[i]] = value
		else:
			LanczosParamsLocal[rangesListText[i]] = value
		
		#Get the potential grid range
		HvParamsLocal["potentialGridSize"] = round(HvParamsLocal["gridSize"] + H2Allowance, 15)
		
		#Get the number of grid points
		HvParamsLocal["gridPoints"] = int(round((HvParamsLocal["gridSize"] / HvParamsLocal["gridSpacing"]) + 1.))
		HvParamsLocal["potentialGridPoints"] = int(round((HvParamsLocal["potentialGridSize"] / HvParamsLocal["potentialGridSpacing"]) + 1.))
		
		#Test if the simulation is the same as the base parameters, if so, don't make a new file for it
		shared_items1 = set(HvParamsLocal.items()) & set(baseHvParamDict.items()) #Taken from http://stackoverflow.com/questions/4527942/comparing-two-dictionaries-in-python
		shared_items2 = set(LanczosParamsLocal.items()) & set(baseLanczosParamDict.items())
		
		#print rangesListText[i], ": ", value
		#print set(HvParamsLocal.items())
		#print set(baseHvParamDict.items())
		#print shared_items1
		#print shared_items1 ^ set(HvParamsLocal.items())
		#print shared_items1 ^ set(baseHvParamDict.items())
		#print len(set(HvParamsLocal.items()))
		#print len(shared_items1)
		#raw_input("Press Enter to continue...")
		
		if (len(shared_items1) == len(set(HvParamsLocal.items()))) and (len(shared_items2) == len(set(LanczosParamsLocal.items()))):
			datafile = open(outDirRange + "/" + "baseParameters.txt", 'w')
			datafile.write("This file indicates that the base parameters are part of this calculation range" + "\n")
			datafile.close()
			print "     %s - identical to value in Base Parameters" % (str(value))
			
		else:
		
			print "     %s" % (str(value))
			
			cnt += 1
			#Get new simulation title
			HvParamsLocal["simulationTitle"] = genSimTitle(simulationTitlePrefix, HvParamsLocal, LanczosParamsLocal)
		
			#Write out data to both range and compilation folders
			HvFileTitle = HvFilePrefix + "_" + HvParamsLocal["simulationTitle"] + ".txt"
			writeHvParams(HvParamsLocal, outDirRange + "/" + HvFileTitle)
			writeHvParams(HvParamsLocal, outDirComp + "/" + HvFileTitle)
			
			LanczosFileTitle = LanzcosFilePrefix + "_" + HvParamsLocal["simulationTitle"] + ".txt"
			writeLanczosParams(LanczosParamsLocal, outDirRange + "/" + LanczosFileTitle)
			writeLanczosParams(LanczosParamsLocal, outDirComp + "/" + LanczosFileTitle)
			
			#Append execution command to shell script
			appendShellFile(executable, nthreads, HvParamsLocal, HvFileTitle, LanczosFileTitle, OutputDirectoryPath, outDirComp + "/" + "Run.sh")
	
print ""
print "%d simulations prepared." % (cnt)
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		

