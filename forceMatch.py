#!/Users/martinmccullagh/anaconda/bin/python


import scipy
import sys
import os
import numpy
import math
import MDAnalysis
import numpy.linalg


# read in command line arguments
inPsfFile = sys.argv[1]
coordDcdFile = sys.argv[2]
forceDcdFile = sys.argv[3]

print "PSF file:", inPsfFile
print "Coord DCD file:", coordDcdFile
print "Force DCD file:", forceDcdFile

# define coordinate max/min and binsize
rMin = 0
rMax = 25
rBinSize = 0.5
nBins = int((rMax-rMin)/rBinSize)
print "Number of bins:", nBins
# average force array
fMagAvg = numpy.zeros(nBins,dtype=float)
# count array
fMagCount = numpy.zeros(nBins,dtype=int)

# start MDAnalysis with a universal
#force universe
force = MDAnalysis.Universe(inPsfFile,forceDcdFile)
#coordinate universe
coord = MDAnalysis.Universe(inPsfFile,coordDcdFile)

#Initial atom selection (select last segment which happens to be our ions)
ionsCoord = coord.atoms.segments[-1]
ionsForce = force.atoms.segments[-1]

#print some general log info
print "Numer of time steps in coordinate trajectory:", len(coord.trajectory)
print "Numer of time steps in force trajectory:", len(force.trajectory)

#initialize distance vector
r = numpy.empty(3,dtype=float)

#iterate over coordinate trajectory
for ts in coord.trajectory:
	# get box dimensions
	box = coord.trajectory.ts.dimensions[:3]
	# compute distance vector in 3D
	for i in range(0,3):
		r[i] = ionsCoord[0].positions[0,i] - ionsCoord[1].positions[0,i]
		# check which image is closest
		if r[i] > box[i]/2:
			r[i] -= box[i]
		elif r[i] < -box[i]/2:
			r[i] += box[i]
	# compute the magnitude of the distance vector (aka the distance)	
	magR = numpy.linalg.norm(r)
	# compute the normal vector in the distance direction
	rHat = r/magR

	# Now work in Force Universe
	naForce = ionsForce[0].positions
	clForce = ionsForce[1].positions
	#compute the projection of the force in the distance direction
	naMagF = numpy.dot(naForce,rHat)	
	clMagF = numpy.dot(clForce,rHat)

	# First compute array index of magR
	magRbin = int((magR-rMin)/rBinSize)
	# add to average force array
	if magRbin<nBins and magRbin>=0:
		fMagAvg[magRbin] += naMagF[0]
		fMagAvg[magRbin] += clMagF[0]
		#add to counts
		fMagCount[magRbin] += 2

	# progress the force trajectory
	if ts.frame < len(force.trajectory):
		force.trajectory.next()	


#complete averaging and pring
for i in range(0,nBins):
	print i*rBinSize,fMagAvg[i]/float(fMagCount[i]), fMagCount[i]

