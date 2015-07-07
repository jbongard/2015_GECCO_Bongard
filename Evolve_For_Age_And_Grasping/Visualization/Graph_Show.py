import pydot
from igraph import *
import igraph

from numpy import *
from pylab import *
from scipy import *
from os.path import *

def PyDot_SaveGraph(graphIndex):

    graph = pydot.Dot(graph_type='graph');

    fileName = '../Code/Modularity/Release/Data/GRN_'+str(graphIndex)+'.dat';

    f = open(fileName,'r');

    line = f.readline();

    line = line.split();
    N = int(line[0]);

    for i in range(0,N):
        node = pydot.Node("gene%d" % i);
        graph.add_node(node);
    
    i = 0;
    for line in f.readlines():
        line = line.split();
        j = 0;
        for element in line:
            element = int(element);
            if ( element != 0 ):
                print i,j,element;
                edge = pydot.Edge("gene%d" % i, "gene%d" % j);
                graph.add_edge(edge);
        
            j = j + 1;
        i = i + 1;
        
    f.close();

    fileName = '../Code/Modularity/Release/Data/GRN_'+str(graphIndex)+'.png';
    graph.write_png(fileName,prog='circo')

def DrawGraph(randIndex,phaseIndex,grnIndex):

    fig = figure(randIndex+1);
    clf();
    
    fileName = '../Code/Modularity/Data/GRN_'+str(randIndex)+'_'+str(phaseIndex)+'_'+str(grnIndex)+'.dat';

    f = open(fileName,'r');

    line = f.readline();

    line = line.split();
    N = int(line[0]);

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

    i = 0;
    for line in f.readlines():
        line = line.split();
        j = 0;
        for element in line:
            element = int(element);
            if ( element != 0 ):
                plot([pos[i,0],pos[j,0]],[pos[i,1],pos[j,1]],'k-');
                if ( i==j ):
                    plot(pos[i,0],pos[i,1],'o',markerfacecolor=[1,1,1],markeredgewidth=3,markersize=18);
        
            j = j + 1;
        i = i + 1;

    for i in range(0,N):
        plot(pos[i,0],pos[i,1],'o',markerfacecolor=[1,1,1],markeredgewidth=1,markersize=18);
        text(pos[i,0],pos[i,1],str(i));
        
    f.close();

    fileName = '../Code/Modularity/Data/GRN_'+str(randIndex)+'_'+str(phaseIndex)+'_'+str(grnIndex)+'.png';
    savefig(fileName);
    
    #show();

def Graph_Modularity(randIndex,graphIndex):

    fileName = '../Code/Modularity/Release/Data/GRN_'+str(randIndex)+'_'+str(graphIndex)+'.dat';

    f = open(fileName,'r');

    line = f.readline();

    line = line.split();
    N = int(line[0]);
    
    g = Graph(N);

    i = 0;
    for line in f.readlines():
        line = line.split();
        j = 0;
        for element in line:
            element = int(element);
            if ( element != 0 ):
                g.add_edges((i,j));
            j = j + 1;
        i = i + 1;
        
    f.close();
    
    q1 = g.community_leading_eigenvector().modularity;
    q2 = g.community_label_propagation().modularity;
    q3 = g.simplify().community_fastgreedy().modularity;

    return( (q1+q2+q3)/3 );

def Graph_Modularity(randSeed,phase,grnIndex,regimeName):

    fileName = '../Code/'+str(regimeName)+'/Data/GRN_'+str(randSeed)+'_'+str(phase)+'_'+str(grnIndex)+'.dat';

    f = open(fileName,'r');

    line = f.readline();

    line = line.split();
    N = int(line[0]);
    
    g = Graph(N);

    i = 0;
    for line in f.readlines():
        line = line.split();
        j = 0;
        for element in line:
            element = int(element);
            if ( element != 0 ):
                g.add_edges((i,j));
            j = j + 1;
        i = i + 1;
        
    f.close();
    
    q1 = g.community_leading_eigenvector().modularity;
    q2 = g.community_label_propagation().modularity;
    q3 = g.simplify().community_fastgreedy().modularity;

    return( (q1+q2+q3)/3 );
    #return(q3);

def Fit_Modularity(runIndex):

    figure(runIndex);
    fileName = '../Code/Modularity/Release/Data/runData_'+str(runIndex)+'.dat';
    f = open(fileName,'r');
    fits = [];
    for line in f.readlines():
        line = line.split();
        for element in line:
            fits.append(float(element));
    f.close();
    plot(fits);

    q = [];
    for i in range(0,len(fits)):
        q.append(Graph_Modularity(runIndex,i));
    plot(q);

def Comparative_Modularity(figNo,firstPhase,lastPhase,regimeName):

    numRuns = 400;

    modValues = zeros((numRuns,2,1),dtype='f');
    modBeforeAfter = zeros((numRuns,2),dtype='f');

    for r in range(0,numRuns):
        #print r,numRuns;
        #for phase in range(0,1+1):
        for grn in range(0,0+1):
            modValues[r,0,grn] = Graph_Modularity(r,firstPhase,grn,regimeName);
            modValues[r,1,grn] = Graph_Modularity(r,lastPhase,grn,regimeName);

    for r in range(0,numRuns):
        x = mean( modValues[r,0,:] );
        y = mean( modValues[r,1,:] );
        modBeforeAfter[r,0] = x;
        modBeforeAfter[r,1] = y;

    figure(figNo);
    
    plot( modBeforeAfter[:,0] , modBeforeAfter[:,1] , 'k.');

    ax = axis();
    minAx = min(ax);
    maxAx = max(ax);
    plot([minAx,maxAx],[minAx,maxAx],'k-');

    meanX = mean( modBeforeAfter[:,0] );
    meanY = mean( modBeforeAfter[:,1] );

    plot(meanX,meanY,'r.');

    stdX = std( modBeforeAfter[:,0] );
    stdY = std( modBeforeAfter[:,1] );

    plot([meanX-stdX,meanX+stdX],[meanY,meanY],'r-');
    plot([meanX,meanX],[meanY-stdY,meanY+stdY],'r-');

    print meanX, meanY;

    mostModularGRN = argmax(modValues[:,1,:]);
    modularityScore = max(modValues[:,1,:]);

    print mostModularGRN, modularityScore;

    title(regimeName);
    
    #DrawGraph(mostModularGRN,1,0);

    fileName = '../Images/Fig'+str(figNo)+'.png';
    savefig(fileName);

def Mean_Connectivity(figNo,phaseIndex,regimeName,numGenes):

    numRuns = 400;

    figure(figNo);
    
    adjMatrix = zeros((numGenes,numGenes),dtype='f');

    for k in range(0,numRuns):

        fileName = '../Code/'+regimeName+'/Data/GRN_'+str(k)+'_'+str(phaseIndex)+'_0.dat';

        f = open(fileName,'r');

        line = f.readline();

        line = line.split();
        N = int(line[0]);

        i = 0;
        for line in f.readlines():
            line = line.split();
            j = 0;
            for element in line:
                element = int(element);
                if ( element != 0 ):
                    adjMatrix[i,j] = adjMatrix[i,j]+1;
                j = j + 1;
            i = i + 1;
            
        f.close();

    print adjMatrix;
    maxVal = adjMatrix.max();
    print maxVal;

    for i in range(0,numGenes):
        for j in range(0,numGenes):
            adjMatrix[i,j] = 1.0 - adjMatrix[i,j]/maxVal;

    print adjMatrix;

    for i in range(0,N-1):
        for j in range(i+1,N):
            adjMatrix[i,j] = adjMatrix[i,j]+adjMatrix[j,i]/2.0;
            adjMatrix[j,i] = 0;

    print adjMatrix;
    
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

    print pos;

    for i in range(0,N-1):
        for j in range(i+1,N):
            cComponent = adjMatrix[i,j];
            if ( cComponent > 1.0 ):
                cComponent = 1.0;
            c = [cComponent,cComponent,cComponent];
            plot([pos[i,0],pos[j,0]],[pos[i,1],pos[j,1]],'k-',color=c);
            print i,j;

    for i in range(0,N):
        plot(pos[i,0],pos[i,1],'ko',markerfacecolor=[1,1,1],markersize=18);

    title(regimeName);
    
    fileName = '../Images/Fig'+str(figNo)+'.png';
    savefig(fileName);
        
# ------------------------------ Main function -------------------------------

#PyDot_SaveGraph(197);

#DrawGraph(1,466);

#Fit_Modularity(1);
#Fit_Modularity(2);

#Comparative_Modularity(3,0,2,'Modularity_3_NoBiasedMutation');
#Comparative_Modularity(2,0,2,'Mod3Phases_2Modules');

#Comparative_Modularity(1,0,1,'Modularity');
#Comparative_Modularity(2,0,2,'Modularity_3');
#Comparative_Modularity(3,0,3,'Mod4Phases_2Modules');
#Comparative_Modularity(4,0,1,'Mod_StartEndDifferent');
#Comparative_Modularity(5,0,1,'Mod_SESimilar');
#Comparative_Modularity(6,0,3,'Mod_Sensor');

#Mean_Connectivity(2,2,'Modularity_3');
#Mean_Connectivity(3,2,'Modularity_3_NoBiasedMutation');

Mean_Connectivity(7,1,'Modularity',10);
Mean_Connectivity(8,2,'Mod3Phases_2Modules',10);
Mean_Connectivity(9,3,'Mod4Phases_2Modules',10);
Mean_Connectivity(10,1,'Mod_StartEndDifferent',10);
Mean_Connectivity(11,1,'Mod_SESimilar',10);
Mean_Connectivity(12,3,'Mod_Sensor',12);

show();
