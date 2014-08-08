#Calculates cross-sections and filtering efficiency from HepMC files
#David Wardrope, 29th July 2014
import sys, re

print "Processing files %s" % (sys.argv[1:])

avgXSec = 0
totalGenEvt = 0
totalFiltEvt = 0
for fName in sys.argv[1:]:
	with open(fName, 'r') as f:
		xsec = 0
		while xsec == 0:
			l = f.readline()
			words = l.split()
			if len(words) > 0:
				if words[0] == 'C':
					xsec = float(words[1])
		print "%s cross-section = %.2f" % (fName, xsec)
		avgXSec += xsec
		f.seek(-1e6, 2)
		l = f.readline()
		nFiltEvt = 0
		while l:
			words = l.split()
			if len(words) > 0:
				if words[0] == 'E':
					nFiltEvt = int(words[1])
			#print l
			l = f.readline()
		print "%s filter efficiency = %d/100000 = %.2f" % (fName, nFiltEvt, nFiltEvt/100000.)
		if nFiltEvt > 0:
			totalFiltEvt += nFiltEvt
			totalGenEvt += 100000

avgXSec /= len(sys.argv[1:])
print "Cross-section = %.4f pb" % (avgXSec)
print "Filter efficiency = %d/%d = %.4f" % (totalFiltEvt, totalGenEvt, totalFiltEvt/(1.*totalGenEvt))
effXSec = (avgXSec*totalFiltEvt)/(1.*totalGenEvt)
print "Effective cross-section = %.4f pb " % (effXSec)
print "Effective integrated luminosity = %.4f fb^{-1}" % (totalFiltEvt/(1000.*effXSec))

