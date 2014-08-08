#Calculates cross-sections and filtering efficiency from HepMC files
#David Wardrope, 29th July 2014
import sys, re

print "Processing files %s" % (sys.argv[1:])

avgXSec = 0
totalGenEvt = 0
for fName in sys.argv[1:]:
	with open(fName, 'r') as f:
		sumEvtWeights = 0.0
		sumNumTrials = 0.0
		for line in f:
			words = line.split()
			if len(words) > 0:
				if words[0] == 'E':
					totalGenEvt += 1
					sumEvtWeights += float(words[13])
					sumNumTrials += float(words[16])
		xsec = sumEvtWeights/sumNumTrials
		print "%s cross-section = %.3f pb" % (fName, xsec)
		avgXSec += xsec

avgXSec /= len(sys.argv[1:])
print "Cross-section = %.4f pb" % (avgXSec)
print "Number of generated events = %d" % (totalGenEvt)
print "Effective integrated luminosity = %.4f fb^{-1}" % (totalGenEvt/(1000.*avgXSec))

