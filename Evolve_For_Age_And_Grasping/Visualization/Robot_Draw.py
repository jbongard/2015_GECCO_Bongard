from numpy import *
from pylab import *
from scipy import *
from os.path import *

NUM_OBJECTS = 4

def Get_Matrix(objectIndex):

    fileName = 'Robot_Matrix_0_'+str(objectIndex)+'.txt';

    f = open(fileName,'r');
    numSegments = 0;
    for line in f.readlines():
        numSegments = numSegments + 1;
    f.close();

    M = zeros((numSegments,5),dtype='f');

    i = 0;
    f = open(fileName,'r');
    for line in f.readlines():
        line = line.split();
        for j in range(0,len(line)):
                print i,j,line[j];
                M[i,j] = float(line[j]);
        i = i + 1;
    f.close();

    return M;

def Draw_Circle(i):

	fig = gcf()

	circle1 = []

	myIndex = 0

	for size in linspace(1.0,0.5,sqrt(NUM_OBJECTS)):

		for position in linspace(2.5,-2.5,sqrt(NUM_OBJECTS)):

			if ( myIndex == i ):

				x = position
				radius = size
				y = radius
				circle1 = Circle((x,y),radius,color='b')

			myIndex = myIndex + 1

	fig.gca().add_artist(circle1)

def Draw_Ground():

	plot([-5,+5],[0,0],'k-');

def Draw_Robot(M):

    ln = size(M,0);
    wd = size(M,1);

    x1=0;
    y1=0;
    x2=0;
    y2=0;

    for i in range(0,ln):
        parentObj = M[i,0];
        length = M[i,1];
        angle = M[i,2];
        x = M[i,3];
        y = M[i,4];
        if ( parentObj == -1 ):
            x1=0;
            y1=0;
        else:
            x1=M[parentObj,3];
            y1=M[parentObj,4];
        x2=x;
        y2=y;
        c = [0,0,0];
	if ( i==3 ):
		c = [0.5,0,0]
	elif ( i==4 ):
		c = [0,0.5,0]
	elif ( i==5 ):
		c = [1,0,0]
	elif ( i==6 ):
		c = [0,1,0]
        plot([x1,x2],[y1,y2],'k-',linewidth=2,color=c);
   
def Robot_Succeeded_In_Environment(o):

	fileName = 'Robot_Matrix_0_' + str(o) + '.txt'

	return ( isfile(fileName) )
 
# Main function ---------------------------------------

figure(1)

for o in range(0,NUM_OBJECTS):

	if ( Robot_Succeeded_In_Environment(o) ):

		subplot(sqrt(NUM_OBJECTS),sqrt(NUM_OBJECTS),o+1)

		M = Get_Matrix(o)
		Draw_Ground()
		Draw_Circle(o)
		Draw_Robot(M)
		ylim([-5,5])
		xlim([-5,5])

show();
