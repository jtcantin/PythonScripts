import sys
import os

#Read in Summary table data
def readData(filename):
	datafile = open(filename, 'r')

	lines = []

	for line in datafile:
		if line.find("  SYM.  ROOT #     EE(eV)      AEL %   %GUESS  #ITER   RESIDUAL    TOTAL ENERGY") != -1:
			lines.append(line)	
			for line2 in datafile:
				if line2.find("Finished calculation") != -1:
					break

				lines.append(line2)


			break
			
	datafile.close()

	if len(lines) == 0:
		sys.exit("ERROR: Data not found!")

	return lines

#Write summary table given list of data lines
def outputData(filename,linelist):
	datafile = open(filename, 'w')
	#Write header
	datafile.write(linelist[0][0])
	datafile.write(linelist[0][1])
	datafile.write(linelist[0][2])

	#Write body
	for filelines in linelist:
		for line in filelines[3:-3]:
			datafile.write(line)

	#Write footer
	datafile.write(linelist[-1][-3])
	datafile.write(linelist[-1][-2])
	datafile.write(linelist[-1][-1])	

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
dirc=sys.argv[1]
topdirclist=os.listdir(dirc)
gooddirc = []

#Gather integer-named directoris
for direct in topdirclist:
	if is_int(direct):
		gooddirc.append(direct)

#Sort gooddirc
gooddirc.sort()

cc3linelist = []
sdt3linelist = []

# Grouping all CC3 summary tables and all sdt3 summary tables into two lists
for direct in gooddirc:
	directPath = "./" + dirc + "/" + direct
	subdirclist = os.listdir(directPath)
	cc3dirclist = []
	sdt3dirclist = []

	#find cc3 output files and sdt3 output files and put them into two respective dirclists
	for subdirc in subdirclist:
		if (subdirc.find("cc3") != -1) and (subdirc.find(".out0") != -1):
			cc3dirclist.append(subdirc)	 
		
		if (subdirc.find("sdt3") != -1) and (subdirc.find(".out0") != -1):
			sdt3dirclist.append(subdirc)
	
	#Sort cc3dirctlist and sdt3dirclist
	cc3dirclist.sort()
	sdt3dirclist.sort()
	
	#print cc3dirclist
	#print sdt3dirclist

	#Read cc3 output summary tables and put in cc3linelist 
	for cc3file in cc3dirclist:
		filename = directPath + "/" + cc3file
		cc3linelist.append(readData(filename))
 
	#Read sdt3 output summary tables and put in sdt3linelist 
	for sdt3file in sdt3dirclist:
		filename = directPath + "/" + sdt3file
		sdt3linelist.append(readData(filename))

#print cc3linelist
#print sdt3linelist

#Process data for CC3 and put all glued lines in summary file
filename = "./" + dirc + "/" + dirc + "_cc3_Summary"
outputData(filename,cc3linelist)


#Process data for sdt3 and put all glued lines in summary file
filename = "./" + dirc + "/" + dirc + "_sdt3_Summary"
outputData(filename,sdt3linelist)

