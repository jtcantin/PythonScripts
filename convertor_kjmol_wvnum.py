import sys
import os

#Read relative eigenvalues
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
		energies.append(line.split()[1])

	return {"ZPE": zpe, "Energies": energies, "Title": title}

#Convert values to wavenumbers from kJ/mol
def kjmolToWvnumData(data):
	energies = data["Energies"]

	h = 6.62606957E-34 #J.s, from NIST 2010
	c = 299792458 #m/s, from NIST 2010
	N_A = 6.02214129E+23 #/mol, from NIST 2010

	wvnumToKjmol = h / 1000.0 * c * 100.0 * N_A
	kjmolTowvnum = 1.0 / wvnumToKjmol

	energies2 = []
	for value in energies:
		energies2.append(float(value) * kjmolTowvnum)

	zpe = data["ZPE"] * kjmolTowvnum
	
	return {"ZPE": zpe, "Energies": energies2, "Title": data["Title"]}



#Write the data to a single file
def outputData(filename,data):
	datafile = open(filename, 'w')
	#Write header
	datafile.write("#The following are the relative energies for various simulations" + "\n")
	datafile.write("#All of the data are in cm-1" + "\n")

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

#Code to test if integer or not.
#Acquired from http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-in-python
def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

#Start of program
dirc = sys.argv[1]
outputFilename = dirc + "/" + sys.argv[2]

#Get the directories and sort them
dirList = [ name for name in os.listdir(dirc) if os.path.isdir(os.path.join(dirc, name)) ] #Taken from http://stackoverflow.com/questions/141291/how-to-list-only-top-level-directories-in-python
dirList.sort()

data = []

for direct in dirList:
	dirPath = "./" + dirc + "/" + direct

	filename = dirPath + "/" + "states_zpe.txt"
	data.append(kjmolToWvnumData(readZPE(filename)))

outputData(outputFilename,data)

