import os;
import time;
import math;
import sys;

def MakeFile(fileIndex,startRun,endRun):

	filename = '../run_'+str(fileIndex)+'.bat'

	f = open(filename,'w')
 
	for i in range(startRun,endRun+1):

		for useEvolvedBodies in range(0,1+1):

			writeString = './Modularity ' + str(100*useEvolvedBodies + i) + ' 1 ' + str(useEvolvedBodies) + '\n';

			f.write(writeString)

	f.close()


MakeFile(1,0,12)
MakeFile(2,13,24)
MakeFile(3,25,36)
MakeFile(4,37,48)
MakeFile(5,49,60)
MakeFile(6,61,72)
MakeFile(7,73,84)
MakeFile(8,85,99)
