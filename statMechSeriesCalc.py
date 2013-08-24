###########################################################################
# statMechSeriesCalc.py
#	This program takes the relative and absolute boundstates as output from the lanczos boundstate calculator,
#  converts the values to wavenumbers, outputs the relative boundstates UP TO the first non-zero error, for both
#  para and ortho hydrogen, then calculates the stat mech properties Z, U, A, S, Cv, Probability of occupation of highest
#  state, and the probability of occupation of each state as a function of temperature, from dT to T_max for ortho, para,
#  true, and normal (3:1 ratio) hydrogen and outputs these to two sets of files for each property. The one set has the 
#  ortho, para, true and normal hydrogen data for each simulation parameter set, while the other file set has the 
#  simulation parameter data for each of ortho, para, true, and normal hydrogen. True hydrogen is that which 
#  incorporates both ortho and para in the partition function with the appropriate spin statistical weighting (3 and 1,
#  respectively) and corresponds to the cooling of room temperature hydrogen in the presence of a ferromagnetic catalyst. 
#  Normal hydrogen assumes a 3:1 ortho:para ratio at all temperatures and corresponds to cooling room temperature
#  hydrogen without a ferromagnetic catalyst that allows ortho-para conversion.
#####################
#	Usage:
#		pydev statMechSeriesCalc.py [orthoFolder] [paraFolder] [outputDirectoryName]
#####################
#
#  First written by Joshua Cantin on 22 August 2013
#
##########################################################################

import sys
import os
import string
import datetime
from sys import argv
import numpy as np
from pylab import *

#Constants from NIST (CODATA 2010)
c = 299792458 #m/s
h = 6.62606957E-34 #J.s
q_e = 1.602176565E-19 #C
N_A = 6.02214129E23 #mol^-1
kB_J = 1.3806488E-23 #J/K
kB = kB_J*N_A/1000 #kJ/mol

#Temperature Range
delta_T = 0.1 #K
T_min = delta_T #K
T_max = 500 #K

#delta_T = sys.argv[4] #K
#T_min = delta_T #K
#T_max = float(sys.argv[5]) #K

###########################################################################
# Functions
###########################################################################

#################################
#Read relative eigenvalues
#################################
def readZPE(filename):
	datafile = open(filename, 'r')

	junk = datafile.readline()

	# Get the simulation title
	line = datafile.readline()
	title = line.split()[9]

	junk = datafile.readline()

	# Get the ground state energy
	line = datafile.readline()
	zpe = float(line.split()[4])

	junk = datafile.readline()

	#Get the energies
	energies = []
	for line in datafile:
		energies.append(float(line.split()[1]))
		
	datafile.close()

	return {"ZPE": zpe, "Energies": energies, "Title": title}

#################################
#Convert values to wavenumbers from kJ/mol
#################################
def kjmolToWvnumData(data):
	energies = data["Energies"]

	h = 6.62606957E-34 #J.s, from NIST 2010
	c = 299792458 #m/s, from NIST 2010
	N_A = 6.02214129E+23 #/mol, from NIST 2010

	wvnumToKjmol = h / 1000.0 * c * 100.0 * N_A
	kjmolTowvnum = 1.0 / wvnumToKjmol

	energies2 = []
	for value in energies:
		energies2.append(value * kjmolTowvnum)

	zpe = data["ZPE"] * kjmolTowvnum
	
	return {"ZPE": zpe, "Energies": energies2, "Title": data["Title"]}

#################################
#Write the data to a single file
#################################
def outputData(filename,data,flag):
	datafile = open(filename, 'w')
	#Write header
	datafile.write("#The following are the relative energies for various simulations" + "\n")

        if flag == "cm":
            datafile.write("#All of the data are in cm-1" + "\n")
        elif flag == "kjmol":
            datafile.write("#All of the data are in kJ/mol" + "\n")
        else:
            print "Error, incorrect flag in outputData: ", flag

	#Write the simulation names
	datafile.write("State")
	for simData in data:
		datafile.write(" " + simData["Title"])

	datafile.write("\n")

	#Write the ZPE
	datafile.write("ZPE")
	for simData in data:
		datafile.write(" " + "%.15E" % simData["ZPE"])

	datafile.write("\n")

	#Determine the largest number of data points
	lengths = []

	for simData in data:
		lengths.append(len(simData["Energies"]))
	lengths.sort()
	maxNum = lengths[-1]

	#Write the state energies
	for i in range(0, maxNum):
		datafile.write("%d" % i)

		for simData in data:
			energies = simData["Energies"]

			if len(energies) <= i:
				datafile.write(" " + "NA")
			else:
				datafile.write(" " + "%.15E" % energies[i])

		datafile.write("\n")

	datafile.close()

#################################
# Function to Claculate Thermodynamic Properties
#################################
#
# !NOTE: the input data should be in cm^-1!
#

def GetH2Thermo(GS_Ref, data): 

	GS_Eng = data["ZPE"]
	print "Original Ground State Energy = ", GS_Eng, " cm^-1"
	GS_Eng -= GS_Ref
	print "New Ground State Energy = ", GS_Eng, " cm^-1"
	print "Adjusted to a Ground State of: ", GS_Ref, " cm^-1"

	relEngWavenumArray = np.array(data["Energies"])

	#print relEngWavenumArray
	#print relEngWavenumArray.size

	EngWavenumArray = relEngWavenumArray+GS_Eng

	#print EngWavenumArray

	EngJouleArray = h*c*100*EngWavenumArray
	EngkJmolArray = EngJouleArray*N_A/1000

	#print EngkJmolArray


	#Do Kronecker product of E*beta
	E_Beta_mat_stor = np.kron(np.transpose(EngkJmolArray), beta_array)

	#Reshape the matrix so that it becomes 2D and not still a 1D array of lists
	E_Beta_mat = np.reshape(E_Beta_mat_stor, (EngkJmolArray.size,-1))
	#print E_Beta_mat
	#print EngkJmolArray.size

	BoltzFactor_mat = np.exp(-1*E_Beta_mat)

	#print BoltzFactor_mat
	#print BoltzFactor_mat[1,:]

	#Calculate the partition function
	PartitonFcn_array = np.sum(BoltzFactor_mat, axis=0)

	#print PartitonFcn_array

	#test = np.array(E_Beta_mat, ndmin=2)
	#print test

	#Calculate the probability of each state
	Prob_mat = np.empty([EngkJmolArray.size, beta_array.size])
	for i in range(0, beta_array.size ):
		Prob_mat[:,i] = BoltzFactor_mat[:,i]/PartitonFcn_array[i]

	#print Prob_mat
	#Calculate the Internal Energy    
	E_Prob_mat = np.empty([EngkJmolArray.size, beta_array.size])
	for i in range(0, EngkJmolArray.size):
		E_Prob_mat[i,:] = Prob_mat[i,:]*EngkJmolArray[i]

	U_array = np.sum(E_Prob_mat, axis=0)

	#print E_Prob_mat
	#print U_array

	#Calculate the average of the energy squared
	E_sq_Prob_mat = np.empty([EngkJmolArray.size, beta_array.size])
	for i in range(0, EngkJmolArray.size):
		E_sq_Prob_mat[i,:] = Prob_mat[i,:]*EngkJmolArray[i]*EngkJmolArray[i]

	E_sq_array = np.sum(E_sq_Prob_mat, axis=0)

	#print E_sq_array

	#Calculate the Variance
	Var_Form_array = E_sq_array+(-1*U_array*U_array)

	#Calculate p*(E-U)^2
	E_U_sq_Prob_mat = np.empty([EngkJmolArray.size, beta_array.size])
	for i in range(0, EngkJmolArray.size):
		E_U_sq_Prob_mat[i,:] = Prob_mat[i,:]*((-1*U_array)+EngkJmolArray[i])**2

	Var_array = np.sum(E_U_sq_Prob_mat, axis=0)

	#print Var_array
	#print Var_Form_array

	#Calulate the Heat Capacity
	Cv_array = Var_array/(kB*T_array*T_array)
	Cv_Form_array = Var_Form_array/(kB*T_array*T_array)

	#print "Cv (T) = Sum(p*(E-U)^2) = ", Cv_array
	#print "Cv (T) = (<E2>-<E>2)/(kB*T^2) = ", Cv_Form_array


	#Caluclate the Free Energy (Helmholtz)
	Helm_array = -1*kB*T_array*np.log(PartitonFcn_array)
	#print "A(T) = -kBT*lnZ = ", Helm_array

	#Calculate the Entropy
	S_array_one = kB*(np.log(PartitonFcn_array) + beta_array*U_array)
	S_array_two = (U_array - Helm_array)/T_array

	#print "S(T) = kB(lnZ + beta*U) = ", S_array_one
	#print "S(T) = (U - A)/T = ", S_array_two

	return {"Z": PartitonFcn_array, "A": Helm_array, "U": U_array, "S": S_array_one, "Var": Var_Form_array,
		"Cv": Cv_Form_array, "E2Av": E_sq_array, "Eng": EngkJmolArray, "EngWvnum" : EngWavenumArray, "ProbMat": Prob_mat}
	
	
def log_fmt(x, pos):
	if x<np.finfo(float).eps:
		ret = 0
	else:
		ret = x*10**(-1*np.floor(np.log10(np.abs(x))))
			
	return "%3.0f" % (ret)

def log_fmt_2(x, pos):
	if x<np.finfo(float).eps:
		ret = 0
	elif (np.floor(np.log10(np.abs(x))))>(-2.1):
		ret = x*10**(-1*(np.floor(np.log10(np.abs(x)))-1))
		#print "here"
	else:
		ret = x*10**(-1*np.floor(np.log10(np.abs(x))))
		#print np.floor(np.log10(np.abs(x)))
			
	return "%3.0f" % (ret)

def R_fmt(x, pos):
	
	return "%3.2f R" % (x)

###########################################################################
# Main Program
###########################################################################

orthoDir = sys.argv[1]
paraDir = sys.argv[2]
outputDir = sys.argv[3] + "_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

#Make output directory
os.system("mkdir ./" + outputDir)

###################################
#Collect and convert ortho data
###################################

dirc = orthoDir


#Get the directories and sort them
dirList = [ name for name in os.listdir(dirc) if os.path.isdir(os.path.join(dirc, name)) ] #Taken from http://stackoverflow.com/questions/141291/how-to-list-only-top-level-directories-in-python
dirList.sort()

#Collect the ortho data and output to a file
orthoDataWvnum = []
orthoDataKjmol = []

for direct in dirList:
	dirPath = "./" + dirc + "/" + direct

	filename = dirPath + "/" + "states_zpe.txt"
	orthoDataKjmol.append(readZPE(filename))
	orthoDataWvnum.append(kjmolToWvnumData(orthoDataKjmol[-1]))

#cm-1 data
outputFilename = outputDir + "/" + "orthoWvnumSummaryZPE.txt"
outputData(outputFilename,orthoDataWvnum,"cm")

#kJ/mol data
outputFilename = outputDir + "/" + "orthoKjmolSummaryZPE.txt"
outputData(outputFilename,orthoDataKjmol,"kjmol")

###################################
#Collect and convert para data
###################################
dirc = paraDir

#Get the directories and sort them
dirList = [ name for name in os.listdir(dirc) if os.path.isdir(os.path.join(dirc, name)) ] #Taken from http://stackoverflow.com/questions/141291/how-to-list-only-top-level-directories-in-python
dirList.sort()

#Collect the para data and output to a file
paraDataWvnum = []
paraDataKjmol = []

for direct in dirList:
	dirPath = "./" + dirc + "/" + direct

	filename = dirPath + "/" + "states_zpe.txt"
	paraDataKjmol.append(readZPE(filename))
	paraDataWvnum.append(kjmolToWvnumData(paraDataKjmol[-1]))

#cm-1 data
outputFilename = outputDir + "/" + "paraWvnumSummaryZPE.txt"
outputData(outputFilename,paraDataWvnum,"cm")

#kJ/mol data
outputFilename = outputDir + "/" + "paraKjmolSummaryZPE.txt"
outputData(outputFilename,paraDataKjmol,"kjmol")

###################################
#Output the simulations to ensure correct ortho and para pairing
###################################
outputFilename = outputDir + "/" + "pairingSummary.txt"

datafile = open(outputFilename, 'w')

datafile.write("#This file lists the ortho and para simulation titles to allow a check that the simulations are properly paired" + "\n")

if len(orthoDataKjmol) != len(paraDataKjmol):
	print "Warning, the number of ortho and para simulations is not equal"

index = 0
for sim in orthoDataKjmol:
	datafile.write("ortho %d:" % index)
	datafile.write(" " + sim["Title"] + "\n")
	
	datafile.write("para  %d:" % index)
	if index < len(paraDataKjmol):
		datafile.write("   " + paraDataKjmol[index]["Title"] + "\n")
	else:
		datafile.write(" " + "NA" + "\n")
	
	datafile.write("---------------------------------------" + "\n")
	
	index += 1

datafile.close()

###################################
#Calculate the thermodynamic properties for each simulation
###################################

T_array_store = np.array(np.arange(T_min, T_max, delta_T))
T_array = np.array([T_min])
T_array = np.append(T_array, T_array_store)

beta_array = pow(kB*T_array, (-1))

numSims = min(len(orthoDataKjmol), len(paraDataKjmol))



#Get thermodynamic data
simOrthoThermoProp = []
simParaThermoProp = []
simTrueThermoProp = []
simNormalThermoProp = []

trueDataWvnum = []
normalDataWvnum = []

for simNum in range(0,numSims):

	print " "
	print "Simulation ", simNum

	#Ortho properties
	#################
	simOrthoThermoProp.append(GetH2Thermo(paraDataWvnum[simNum]["ZPE"], orthoDataWvnum[simNum]))
	
	#Para properties
	#################
	simParaThermoProp.append(GetH2Thermo(paraDataWvnum[simNum]["ZPE"], paraDataWvnum[simNum]))
	
	#Generate True Hydrogen energy levels
	#################
	trueDataWvnum.append({"ZPE": paraDataWvnum[simNum]["ZPE"], "Energies": [], 
		"Title": paraDataWvnum[simNum]["Title"]})
		
	#Add para energies
	trueDataWvnum[-1]["Energies"] += list(simParaThermoProp[-1]["EngWvnum"])
	
	#Add ortho energies three times
	trueDataWvnum[-1]["Energies"] += list(simOrthoThermoProp[-1]["EngWvnum"])
	trueDataWvnum[-1]["Energies"] += list(simOrthoThermoProp[-1]["EngWvnum"])
	trueDataWvnum[-1]["Energies"] += list(simOrthoThermoProp[-1]["EngWvnum"])
	
	
	#Sort the energies
	trueDataWvnum[-1]["Energies"].sort()
	
	#print trueDataWvnum[-1]["Energies"][7:10]
	
	#True properties
	#################
	simTrueThermoProp.append(GetH2Thermo(paraDataWvnum[simNum]["ZPE"], trueDataWvnum[simNum]))
	
	#Normal properties -> NOT DONE YET!
	#################
	

###################################
#Plot data
###################################

prob_limit = 20

for simNum in range(0,numSims):

	#Plot the Heat Capacity relative to R
	figure(simNum+1, figsize=(9,7))
	ax = subplot(111)
	plot(T_array, simParaThermoProp[simNum]["Cv"]/kB, label="p-H$_{2}$")
	plot(T_array, simOrthoThermoProp[simNum]["Cv"]/kB, label="o-H$_{2}$")
	plot(T_array, simTrueThermoProp[simNum]["Cv"]/kB, label="True H$_{2}$")
	#plot(T_array, spinless_prop["Cv"]/kB, label="Spinless H$_{2}$", color="orange")
	legend(loc='upper left')
	xlabel("T (K)", fontsize=17)
	ylabel("C$_V$ ", fontsize=17)
	title("Heat Capacity"+"\n"+"Of One Hydrogen Molecule in Small Cage 11 with a"+"\n"+"Quadrupole Pairwise Potential", fontsize=20)
	ax.xaxis.set_minor_locator(MultipleLocator(1.25))
	ax.yaxis.set_major_locator(MultipleLocator(0.25))
	tick_params(labelsize=14, length=5, which='major')
	tick_params(labelsize=14, length=3, which='minor', direction='in')
	ax.yaxis.set_major_formatter(FuncFormatter(R_fmt))
	#ticklabel_format(style="sci", scilimits=(-1,3))
	xlim(0,prob_limit)
	ylim(0,1)
	tight_layout()
	


#Plot the Heat Capacity relative to R for para
figure(numSims+1, figsize=(9,7))
ax = subplot(111)
for simNum in range(0,numSims):
	plot(T_array, simParaThermoProp[simNum]["Cv"]/kB, label="%d Theta/Phi Points" % (simNum*10 + 10))
#plot(T_array, spinless_prop["Cv"]/kB, label="Spinless H$_{2}$", color="orange")
legend(loc='upper left')
xlabel("T (K)", fontsize=17)
ylabel("C$_V$ ", fontsize=17)
title("Heat Capacity"+"\n"+"Of para-Hydrogen in Small Cage 11 "+"\n"+"Varying the Number of Theta/Phi Points", fontsize=20)
ax.xaxis.set_minor_locator(MultipleLocator(1.25))
ax.yaxis.set_major_locator(MultipleLocator(0.1))
tick_params(labelsize=14, length=5, which='major')
tick_params(labelsize=14, length=3, which='minor', direction='in')
ax.yaxis.set_major_formatter(FuncFormatter(R_fmt))
#ticklabel_format(style="sci", scilimits=(-1,3))
xlim(0,prob_limit)
ylim(0,0.4)
tight_layout()

#Plot the Heat Capacity relative to R for ortho
figure(numSims+2, figsize=(9,7))
ax = subplot(111)
for simNum in range(0,numSims):
	plot(T_array, simOrthoThermoProp[simNum]["Cv"]/kB, label= orthoDataWvnum[simNum]["Title"][44:46] + " Theta/Phi Points")
#plot(T_array, spinless_prop["Cv"]/kB, label="Spinless H$_{2}$", color="orange")
legend(loc='upper left')
xlabel("T (K)", fontsize=17)
ylabel("C$_V$ ", fontsize=17)
title("Heat Capacity"+"\n"+"Of ortho-Hydrogen in Small Cage 11 "+"\n"+"Varying the Number of Theta/Phi Points", fontsize=20)
ax.xaxis.set_minor_locator(MultipleLocator(1.25))
#ax.yaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_major_locator(MultipleLocator(0.025))
tick_params(labelsize=14, length=5, which='major')
tick_params(labelsize=14, length=3, which='minor', direction='in')
ax.yaxis.set_major_formatter(FuncFormatter(R_fmt))
#ticklabel_format(style="sci", scilimits=(-1,3))
xlim(0,prob_limit)
ylim(0,1)
#xlim(6,8)
#ylim(0.62,0.67)
tight_layout()

#Plot the Heat Capacity relative to R for ortho for T=7.00
figure(numSims+4, figsize=(9,7))
ax = subplot(111)
x_array = []
y_array = []
for simNum in range(0,numSims):
	x_array.append(float(orthoDataWvnum[simNum]["Title"][44:46]))
	y_array.append(simOrthoThermoProp[simNum]["Cv"][69]/kB)
	#plot(T_array, simOrthoThermoProp[simNum]["Cv"]/kB, label= orthoDataWvnum[simNum]["Title"][44:46] + " ")
#plot(T_array, spinless_prop["Cv"]/kB, label="Spinless H$_{2}$", color="orange")
plot(x_array, y_array)
legend(loc='upper left')
xlabel("Theta/Phi Points", fontsize=17)
ylabel("C$_V$ ", fontsize=17)
title("Heat Capacity"+"\n"+"Of ortho-Hydrogen in Small Cage 11  at 7K"+"\n"+"Varying the Number of Theta/Phi Points", fontsize=20)
ax.xaxis.set_minor_locator(MultipleLocator(10))
#ax.yaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_major_locator(MultipleLocator(0.01))
tick_params(labelsize=14, length=5, which='major')
tick_params(labelsize=14, length=3, which='minor', direction='in')
ax.yaxis.set_major_formatter(FuncFormatter(R_fmt))
#ticklabel_format(style="sci", scilimits=(-1,3))
xlim(0,80)
ylim(0.63,0.67)
#xlim(6,8)
#ylim(0.62,0.67)
tight_layout()

#Plot the Heat Capacity relative to R for true hydrogen
figure(numSims+3, figsize=(9,7))
ax = subplot(111)
for simNum in range(0,numSims):
	plot(T_array, simTrueThermoProp[simNum]["Cv"]/kB, label="%d Theta/Phi Points" % (simNum*10 + 10))
#plot(T_array, spinless_prop["Cv"]/kB, label="Spinless H$_{2}$", color="orange")
legend(loc='upper left')
xlabel("T (K)", fontsize=17)
ylabel("C$_V$ ", fontsize=17)
title("Heat Capacity"+"\n"+"Of True Hydrogen in Small Cage 11 "+"\n"+"Varying the Number of Theta/Phi Points", fontsize=20)
ax.xaxis.set_minor_locator(MultipleLocator(1.25))
ax.yaxis.set_major_locator(MultipleLocator(0.125))
tick_params(labelsize=14, length=5, which='major')
tick_params(labelsize=14, length=3, which='minor', direction='in')
ax.yaxis.set_major_formatter(FuncFormatter(R_fmt))
#ticklabel_format(style="sci", scilimits=(-1,3))
xlim(0,prob_limit)
ylim(0,0.5)
tight_layout()

	

show()




