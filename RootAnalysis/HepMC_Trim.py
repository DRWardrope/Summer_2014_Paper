#Trims HepMC files to an event range we want
#Yik Tung Ho, 27th November 2014
#Modified from cross section calculation script by
#David Wardrope, 29th July 2014
import sys, re

#print "Processing files %s" % (sys.argv[1:])

totalGenEvt = 0
startRange = False
for fName in sys.argv[1:]:
	with open(fName, 'r') as f:
		for line in f:
			words = line.split()
			if len(words) > 0:
				if words[0] == 'HepMC::Version':
					print line,
				if words[0] == 'HepMC::IO_GenEvent-START_EVENT_LISTING':
					print line,
				if words[0] == 'HepMC::IO_GenEvent-END_EVENT_LISTING':
					print line,
				if words[0] == 'E':
					totalGenEvt += 1
				if totalGenEvt == 0:
					startRange = True
				if startRange:	
					if totalGenEvt == 101:
						#startRange = False
						print 'HepMC::IO_GenEvent-END_EVENT_LISTING',
						sys.exit()
					else:
						print line,
