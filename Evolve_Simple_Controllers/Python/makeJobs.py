import os;
import time;
import math;
import sys;

f = [ [] , [] ]

filename = '../runBodies0.bat'
f[0] = open(filename,'w')

filename = '../runBodies1.bat'
f[1] = open(filename,'w')

for i in range(300,399+1):

	filename = '../runBodies' + str(i%2) + '.bat';

	writeString = './Modularity ' + str(i) + '\n';

	f[i%2].write(writeString)

f[0].close()
f[1].close()
