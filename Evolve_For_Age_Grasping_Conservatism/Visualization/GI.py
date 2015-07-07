from numpy import *
from pylab import *
from scipy import *
from os.path import *

numRuns = 500;

def Get_Matrix():

    fileName = 'Robot_Matrix.txt';

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

def Draw_Robot(M,k,numSteps,markerSize):

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
        c = (1.0 - (double(k)/double(numSteps)))*0.9;
            
        plot([x1,x2],[y1,y2],'k-',linewidth=1*(numSteps-k)+1,color=[c,c,c]);

        if ( (i>=(ln-2)) & (k==(numSteps-1)) ):
            plot([x2,x2],[y2,y2],'ko',markersize=markerSize,markerfacecolor=[0,0,0]);

def FileExists(fileName):

    return( exists(fileName) );

def GetRobotMatrix(fileName):

    i = 7;
    j = 5;
    k = 0;
    
    f = open(fileName,'r');
    for line in f.readlines():
        if ( len(line)==1 ):
            k = k + 1;
    f.close();

    M = zeros((i,j,k),dtype='f');

    i = 0;
    k = 0;
    f = open(fileName,'r');
    for line in f.readlines():
        if ( len(line)==1 ):
            i = 0;
            k = k + 1;
        else:
            line = line.split();
            for j in range(0,len(line)):
                M[i,j,k] = float(line[j]);
            i = i + 1;
            
    f.close();
    
    return(M);

def Draw_Graph(N,adjMatrix):

    pos = zeros((N,2),dtype='f');
    stepSize = 2*3.14159/N;
    i = 0;
    while ( i<N ):
        angle = i*stepSize;
        x = 1.0*sin(angle);
        y = 1.0*cos(angle);
        pos[i,0] = x;
        pos[i,1] = y;
        i = i + 1;

    for i in range(0,5+1):
        for j in range(6,13+1):
            cComponent = 1-adjMatrix[i,j];
            c = [cComponent,cComponent,cComponent];
            if ( cComponent<1.0 ):
                plot([pos[i,0],pos[j,0]],[pos[i,1],pos[j,1]],'k-',color=c);

    for i in range(0,4+1):
        for j in range(i+1,5+1):
            cComponent = 1-adjMatrix[i,j];
            c = [cComponent,cComponent,cComponent];
            if ( cComponent<1.0 ):
                plot([pos[i,0],pos[j,0]],[pos[i,1],pos[j,1]],'k-',color=c);

    for i in range(6,12+1):
        for j in range(i+1,13+1):
            cComponent = 1-adjMatrix[i,j];
            c = [cComponent,cComponent,cComponent];
            if ( cComponent<1.0 ):
                plot([pos[i,0],pos[j,0]],[pos[i,1],pos[j,1]],'k-',color=c);
                
    for i in range(0,N):
        if ( i<=5 ):
            plot(pos[i,0],pos[i,1],'ko',markerfacecolor=[1,1,1],markersize=18);
        else:
            plot(pos[i,0],pos[i,1],'ko',markerfacecolor=[0.5,0.5,0.5],markersize=18);

    xticks([],[]);
    yticks([],[]);
    
def Read_Network(N,adjMatrix,runIndex,panelIndex):

    fileName='../Code/Mod_Shaping_AB/Data/GRN_'+str(runIndex)+'_'+str(panelIndex)+'.dat';

    f = open(fileName,'r');

    line = f.readline();

    localAdjMatrix = zeros((N,N),dtype='f');
    
    for line in f.readlines():
        line = line.split();
        fromNeuron = int(line[0]);
        toNeuron = int(line[1]);
        if ( (fromNeuron>-1) & (toNeuron>-1) & (localAdjMatrix[fromNeuron,toNeuron]==0) ):
            localAdjMatrix[fromNeuron,toNeuron] = 1.0;
        
    f.close();
    
    return adjMatrix+localAdjMatrix;

def Get_RunData(fileName,numFinishers,gens,mods,i):

    f = open(fileName,'r');
    j = 0;
    for line in f.readlines():
        line = line.split();
        gens[numFinishers[j],j] = double(line[2]);
        armArm = double(line[5]);
        handHand = double(line[6]);
        armHand = double(line[7]);
        if ( armHand == 0.0 ):
            armHand = 0.01;
        mods[numFinishers[j],j] = ((armArm+handHand)/2.0)/armHand;
        
        numFinishers[j] = numFinishers[j] + 1.0;
        j = j + 1;
        
    f.close();
    
    return( numFinishers,gens,mods );

def Add_Line(v,lw):

    xSpacing = linspace(1,len(v),len(v));    
    plot(xSpacing,v,color=[0,0,0],linewidth=lw);

def Add_Vertical_Lines(x1,x2,xSpacing,y1,y2,c,lw):

    for j in range(x1,x2,xSpacing):
        plot([j,j],[y1,y2],color=[c,c,c],linewidth=lw);

def Add_Lines(numConditions,fn,figIndex,lineWidth):

    numFinishers = zeros(numConditions,dtype='f');
    gens = zeros((numRuns,numConditions),dtype='f');
    meanGens = zeros(numConditions,dtype='f');
    serrGens = zeros(numConditions,dtype='f');
    mods = zeros((numRuns,numConditions),dtype='f');
    meanMods = zeros(numConditions,dtype='f');
    serrMods = zeros(numConditions,dtype='f');
        
    for i in range(0,numRuns):
        fileName = '../Code/'+fn+'/Data/runData_'+str(i)+'.dat';
        [numFinishers,gens,mods] = Get_RunData(fileName,numFinishers,gens,mods,i);
    
    for j in range(0,numConditions):
        meanGens[j] = gens[0:numFinishers[j],j].mean();
        #meanGens[j] = gens[0:numRuns,j].mean();
        serrGens[j] = gens[0:numFinishers[j],j].std()/sqrt(numFinishers[j]);
        meanMods[j] = mods[0:numFinishers[j],j].mean();
        serrMods[j] = mods[0:numFinishers[j],j].std()/sqrt(numFinishers[j]);
    
    figure(figIndex);
    Add_Line(numFinishers,lineWidth);
    
    figure(figIndex+1);
    Add_Line(meanGens,lineWidth);
    #Add_Line(meanGens+serrGens,1);
    #Add_Line(meanGens-serrGens,1);
    
    figure(figIndex+2);
    Add_Line(meanMods,lineWidth);   
    Add_Line(meanMods+serrMods,1);
    Add_Line(meanMods-serrMods,1);
    
def Fig1(figIndex):

    figure(figIndex,figsize=(4,10));
    #figure(figIndex);
    
    k = 1;

    for j in range(0,15+1):
                
        for i in range(0,3+1):

            subplot(16,4,k);
            ax=gca();
                    
            xticks([],[]);
            yticks([],[]);

            if ( i==0 ):
                labelStr = '';
                if ( j<15 ):
                    labelStr = labelStr + 'IC' + str(j+1);
                else:
                    labelStr = labelStr + 'Final';
                ylabel(labelStr,fontsize=16);

            if ( j==15 ):
                labelStr = 'Env'+str(i+1);
                xlabel(labelStr,fontsize=16);

            plot([-4.6,4.6],[0,0],color=[0.5,0.5,0.5]);

            if ( i==0 ):
                circ=Circle((-3.268,-0.853),radius=1.207,facecolor=[1,1,1]);
            elif ( i==1 ):
                circ=Circle((+2.793,-0.379),radius=0.536,facecolor=[1,1,1]);
            elif ( i==2 ):
                circ=Circle((+3.268,-0.853),radius=1.207,facecolor=[1,1,1]);
            elif ( i==3 ):
                circ=Circle((-2.793,-0.379),radius=0.536,facecolor=[1,1,1]);  
            ax.add_patch(circ);

            fileName = '../Code/Mod_IC/Data/robot_'+str(i)+'_'+str(j)+'.dat';
            M = GetRobotMatrix(fileName);
            #print M[:,:,0];
            Draw_Robot(M[:,:,0],1,2,2);
            #raw_input('');
            
            axis([-4.6,4.6,-2.1,4.1]);
        
            k = k + 1;

            subplots_adjust(left=0.11,bottom=0.03,right=0.97,top=0.99);

    savefig('../Images/Fig'+str(figIndex)+'.png');
    savefig('../Images/Fig'+str(figIndex)+'.eps');
    
def Fig2(figIndex):

    runIndex = 13;

    figure(figIndex);
            
    for l in range(0,3+1):
        fileName = '../Code/Mod_Shaping_AB/Data/robot_'+str(runIndex)+'_'+str(l)+'.dat';
        M = GetRobotMatrix(fileName);

        subplot(2,2,l+1);
        ax=gca();
        
        plot([-4.6,4.6],[0,0],color=[0.5,0.5,0.5]);
        
        numSteps = size(M,2);

        if ( l==0 ):
            title('(a) Environment 1');
            circ=Circle((-3.268,-0.853),radius=1.207,facecolor=[1,1,1]);
        elif ( l==1 ):
            title('(b) Environment 2');
            circ=Circle((+2.793,-0.379),radius=0.536,facecolor=[1,1,1]);
        elif ( l==2 ):
            title('(c) Environment 3');
            circ=Circle((+3.268,-0.853),radius=1.207,facecolor=[1,1,1]);
        elif ( l==3 ):
            title('(d) Environment 4');
            circ=Circle((-2.793,-0.379),radius=0.536,facecolor=[1,1,1]);  
        ax.add_patch(circ);
                
        for k in range(0,numSteps):

            Draw_Robot(M[:,:,k],k,numSteps,5);

        xticks([],[]);
        yticks([],[]);
           
        axis([-4.6,4.6,-2.1,4.1]);

    subplots_adjust(left=0.01,bottom=0.02,right=0.98,top=0.94);

    savefig('../Images/Fig'+str(figIndex)+'.png');
    savefig('../Images/Fig'+str(figIndex)+'.eps');

def Fig3(figIndex):

    runIndex = 13;
        
    N = (3+2*2)*2;

    figure(figIndex);

    adjMatrix = zeros((N,N,4),dtype='f');

    for i in range(0,3+1):

        subplot(2,2,i+1);
        
        adjMatrix[:,:,i] = Read_Network(N,adjMatrix[:,:,i],runIndex,i);

        for k in range(0,N):
            for l in range(0,N):
                adjMatrix[k,l,i] = (adjMatrix[k,l,i] + adjMatrix[l,k,i])/2.0;
                adjMatrix[l,k,i] = (adjMatrix[l,k,i] + adjMatrix[k,l,i])/2.0;
                
        Draw_Graph(N,adjMatrix[:,:,i]);

        if ( i==0 ):
            title('(a) Environment 1');
        elif ( i==1 ):
            title('(b) Environment 2');
        elif ( i==2 ):
            title('(c) Environment 3');
        elif ( i==3 ):
            title('(d) Environment 4');

    subplots_adjust(left=0.01,bottom=0.02,right=0.98,top=0.94);

    savefig('../Images/Fig'+str(figIndex)+'.png');
    savefig('../Images/Fig'+str(figIndex)+'.eps');
    
def Fig4(figIndex):

    N = (3+2*2)*2;

    figure(figIndex);

    adjMatrix = zeros((N,N,4),dtype='f');

    for i in range(0,3+1):

        for j in range(0,numRuns):
            adjMatrix[:,:,i] = Read_Network(N,adjMatrix[:,:,i],j,i);

        for k in range(0,N):
            for l in range(0,N):
                adjMatrix[k,l,i] = adjMatrix[k,l,i] + adjMatrix[l,k,i];
                adjMatrix[l,k,i] = adjMatrix[l,k,i] + adjMatrix[k,l,i];
        
    adjMatrix = adjMatrix - adjMatrix.min();
    adjMatrix = adjMatrix/adjMatrix.max();
    
    for i in range(0,3+1):

        subplot(2,2,i+1);
    
        Draw_Graph(N,adjMatrix[:,:,i]);

        if ( i==0 ):
            title('(a) Environment 1');
        elif ( i==1 ):
            title('(b) Environment 2');
        elif ( i==2 ):
            title('(c) Environment 3');
        elif ( i==3 ):
            title('(d) Environment 4');

    subplots_adjust(left=0.01,bottom=0.02,right=0.98,top=0.94);

    savefig('../Images/Fig'+str(figIndex)+'.png');
    savefig('../Images/Fig'+str(figIndex)+'.eps');
    
def Fig5_6(figIndex,figFileIndex,expRegime,titleStr,xStep,xLabel):

    numEnvs = 4;
    numICs = 15;
    numConditions = numEnvs * numICs;

    figure(figIndex);
    Add_Vertical_Lines(0,numConditions,1,0,numRuns+10,0.9,1);
    Add_Vertical_Lines(0,numConditions,xStep,0,numRuns+10,0.8,2);
    
    figure(figIndex+1);
    Add_Vertical_Lines(0,numConditions,1,0,1200000,0.9,1);
    Add_Vertical_Lines(0,numConditions,xStep,0,1200000,0.8,2);

    figure(figIndex+2);
    Add_Vertical_Lines(0,numConditions,1,1,2.8,0.9,1);
    Add_Vertical_Lines(0,numConditions,xStep,1,2.8,0.8,2);
    
    Add_Lines(numConditions,expRegime,figIndex,6);
    
    Add_Lines(numConditions,expRegime+'_AnyAtt',figIndex,4);
    
    Add_Lines(numConditions,expRegime+'_NoAtt',figIndex,2);
    
    if ( figIndex==5 ):
        xTicks = linspace(numEnvs,numConditions,numICs);
        xLabels = zeros(numICs,dtype='i');
        for i in range(0,numICs):
            xLabels[i] = i+1;
    else:
        xTicks = linspace(numICs,numConditions,numEnvs);
        xLabels = zeros(numEnvs,dtype='i');
        for i in range(0,numEnvs):
            xLabels[i] = i+1;
            
    figure(figIndex);
    axis([0,numConditions,0,numRuns+10]);
    xticks(xTicks,xLabels,fontsize=18);
    xlabel(xLabel,fontsize=18);
    yticks(fontsize=18);
    ylabel('Successful Runs',fontsize=18);
    title('(a) '+titleStr,fontsize=18);
    subplots_adjust(left=0.14,bottom=0.10,right=0.97,top=0.94);
    savefig('../Images/Fig'+str(figFileIndex)+'a.png');
    savefig('../Images/Fig'+str(figFileIndex)+'a.eps');
    
    figure(figIndex+1);
    axis([0,numConditions-1,0,1200000]);
    xticks(xTicks,xLabels,fontsize=18);
    xlabel(xLabel,fontsize=18);
    yTicks = linspace(0,1200000,7);
    yLabels = ['0','20K','40K','60K','80K','100K','120K'];
    yticks(yTicks,yLabels,fontsize=18);
    ylabel('Elapsed Generations',fontsize=18);
    title('(b) '+titleStr,fontsize=18);
    subplots_adjust(left=0.14,bottom=0.10,right=0.97,top=0.94);
    savefig('../Images/Fig'+str(figFileIndex)+'b.png');
    savefig('../Images/Fig'+str(figFileIndex)+'b.eps');
    
    figure(figIndex+2);
    axis([0,numConditions-1,1,2.8]);
    xticks(xTicks,xLabels,fontsize=18);
    xlabel(xLabel,fontsize=18);
    yticks(fontsize=18);
    ylabel('Modularity',fontsize=18);
    title('(c) '+titleStr,fontsize=18);
    subplots_adjust(left=0.14,bottom=0.10,right=0.97,top=0.94);
    savefig('../Images/Fig'+str(figFileIndex)+'c.png');
    savefig('../Images/Fig'+str(figFileIndex)+'c.eps');
    
# Main function ---------------------------------------

numEnvs = 4;
numICs = 15;

#Fig1(1);
#Fig2(2);
#Fig3(3);
#Fig4(4);

#titleStr = 'Add Attractors, then Widen Basin of Attraction';
#xLabel = 'Initial Condition';
#Fig5_6(5,5,'Mod_Shaping_AB',titleStr,numEnvs,xLabel);

titleStr = 'Widen Basin of Attraction, Then Add Attractors';
xLabel = 'Environment';
Fig5_6(8,6,'Mod_Shaping_BA',titleStr,numICs,xLabel);

show();
