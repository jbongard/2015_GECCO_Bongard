import os;
import time;
import math;
import sys;
import scipy;
import scipy.stats;

def CompareStatistics(startSeed1,startSeed2,index):

	control = []
	experimental = []

        for i in range(0,99+1):

        	randomSeed = startSeed1 + i 

               	fileName = '../Data/results_' + str(randomSeed) + '.txt';

		if ( os.path.isfile(fileName) ):

			f = open(fileName)

			l = f.readline()

			l = l.split()

			# print l

			control.append( float(l[index]) )

			f.close()	

                randomSeed = startSeed2 + i

                fileName = '../Data/results_' + str(randomSeed) + '.txt'; 

                if ( os.path.isfile(fileName) ):

                        f = open(fileName)

                        l = f.readline()

                        l = l.split()

                        # print l 

                        experimental.append( float(l[index]) )

                        f.close()

	print sum(control)/len(control), len(control)

	print sum(experimental)/len(experimental), len(experimental)

	print scipy.stats.mannwhitneyu(control,experimental)

def CompareStats(startSeed1,startSeed2):

	print 'comparing run sets '+str(startSeed1)+' and '+str(startSeed2)+'.'

	print 'comparision of modularity:'
	CompareStatistics(startSeed1,startSeed2,1)
	print ''

	# print 'comparision of grasping ability:' 
	# CompareStatistics(startSeed1,startSeed2,3)
	# print ''

        print 'comparision of attractorAttainment:'
        CompareStatistics(startSeed1,startSeed2,5)
        print ''

	print 'comparision of angle similarity:'
	CompareStatistics(startSeed1,startSeed2,7)
	print ''

        print 'comparision of node differences:'
        CompareStatistics(startSeed1,startSeed2,9)
        print ''

        print 'comparision of pose quality:'
        CompareStatistics(startSeed1,startSeed2,11)
        print ''

# ------------ Main function -------------------

CompareStats(0,100)

