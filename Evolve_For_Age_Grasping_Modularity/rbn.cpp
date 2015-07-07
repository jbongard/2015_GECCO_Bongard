#ifndef _RANDOM_BOOLEAN_NETWORK_CPP
#define _RANDOM_BOOLEAN_NETWORK_CPP

#include "rbn.h"
#include "math.h"

extern int		NUM_NODES;
extern int		START_NUM_CONNECTIONS;
extern int		NUM_UPDATES;
extern int		NUM_NODES_PER_BODY_PART;
extern double		MUTATION_PROBABILITY;

RANDOM_BOOLEAN_NETWORK::RANDOM_BOOLEAN_NETWORK(void) {

	numNodes = NUM_NODES;

	modularity = 1.0;

	EdgeList_Initialize();

	NodeValues_Initialize();
}

RANDOM_BOOLEAN_NETWORK::~RANDOM_BOOLEAN_NETWORK(void) {

	if ( nodeValues ) {

		delete nodeValues;
		nodeValues = NULL;
	}

	if ( edgeList ) {

		delete edgeList;
		edgeList = NULL;
	}
}

void RANDOM_BOOLEAN_NETWORK::CopyFrom( RANDOM_BOOLEAN_NETWORK *parent ) {

	edgeList->CopyFrom( parent->edgeList );

	numEdges = parent->numEdges;

	modularity = parent->modularity;
}

double RANDOM_BOOLEAN_NETWORK::Distance_From_Target_Pattern(int objectIndex) {

	double distance = 0.0;

	int i = NUM_UPDATES-1;

	for (int j=0 ; j < NUM_NODES ; j++ ) {

		if ( (objectIndex==0) || (objectIndex==1) ) {

			for ( int j=0; j < int(double(NUM_NODES)/2.0) ; j++ )

				distance = distance + fabs( nodeValues->Get(i,j) - -1 );
		}
		else {

                        for ( int j=0; j < int(double(NUM_NODES)/2.0) ; j++ )

                                distance = distance + fabs( nodeValues->Get(i,j) - +1 );
		}

                if ( (objectIndex==0) || (objectIndex==2) ) {

                        for ( int j=int(double(NUM_NODES)/2.0) ; j < NUM_NODES ; j++ )

                                distance = distance + fabs( nodeValues->Get(i,j) - -1 );
                }
                else {

                        for ( int j=int(double(NUM_NODES)/2.0) ; j < NUM_NODES ; j++ )

                                distance = distance + fabs( nodeValues->Get(i,j) - +1 );
                }

	}

	return distance;
}

double RANDOM_BOOLEAN_NETWORK::Fraction_Of_Nodes_Unchanging(void) {

	int numNodesUnchanging = 0;

	for (int j=0; j < NUM_NODES; j++ )

		if ( nodeValues->Get(NUM_UPDATES-1,j) == nodeValues->Get(NUM_UPDATES-2,j) )

			numNodesUnchanging++;

	return double(numNodesUnchanging) / double(NUM_NODES);

}

int RANDOM_BOOLEAN_NETWORK::Get_Angle_Index(int bodyPartIndex) {

	if ( NUM_NODES_PER_BODY_PART == 2 )

		return Get_Angle_Index_2(bodyPartIndex);

	else // NUM_NODES_PER_BODY_PART == 4

		return Get_Angle_Index_4(bodyPartIndex);
}

double RANDOM_BOOLEAN_NETWORK::Get_Modularity(void) {

	return modularity;
}

int  RANDOM_BOOLEAN_NETWORK::Get_NodeValue(int nodeIndex) {

	return nodeValues->Get(NUM_UPDATES-1,nodeIndex);
}

int  RANDOM_BOOLEAN_NETWORK::In_Point_Attractor(void) {

	int inPointAttractor = true;

	int j=0;

	while ( (j<NUM_NODES) && (inPointAttractor==true) ) {

		if ( nodeValues->Get(NUM_UPDATES-1,j) != nodeValues->Get(NUM_UPDATES-2,j) )

			inPointAttractor = false;
		else
			j++;
	}

	return inPointAttractor;
}

void RANDOM_BOOLEAN_NETWORK::Mutate(void) {

        int mutated = false;

        for (int u=0;u<numNodes;u++)

                if ( (float(rand())/float(RAND_MAX)) < MUTATION_PROBABILITY ) {

                        MutateGene(u);
                        mutated = true;
                }

	Compute_Modularity();
}

void RANDOM_BOOLEAN_NETWORK::Perturb_Node(int nodeIndex) {

	if ( nodeValues->Get(0,nodeIndex) == -1 )

		nodeValues->Set(0,nodeIndex,+1);
	else
		nodeValues->Set(0,nodeIndex,-1);
}

void RANDOM_BOOLEAN_NETWORK::Print(void) {

        int i = NUM_UPDATES-1;

        for (int j=0; j < NUM_NODES ; j++ ) {

		if ( j==( NUM_NODES_PER_BODY_PART * 3 ) )

			printf(" ");

                if ( nodeValues->Get(i,j) > 0 )

                        printf("P");
                else
                        printf("N");
	}

	printf("\n");

}

void RANDOM_BOOLEAN_NETWORK::Set_Initial_Conditions(int environmentIndex) {

	for (int n=0 ; n<NUM_NODES ; n=n+2 ) { // Set size sensors

		if ( (environmentIndex==0) || (environmentIndex==1) )

			nodeValues->Set(0,n,+1); // object is large

		else
			nodeValues->Set(0,n,-1); // object is small
	}

        for (int n=1 ; n<NUM_NODES ; n=n+2 ) { // Set position sensors

                if ( (environmentIndex==0) || (environmentIndex==2) ) 

                        nodeValues->Set(0,n,+1); // object is on the right
                
                else 
                        nodeValues->Set(0,n,-1); // object is on the left
        }
}

void RANDOM_BOOLEAN_NETWORK::Update(void) {

        int timeStep = 1;

        while ( timeStep<NUM_UPDATES ) {

                EdgeList_RegulateGenes(timeStep);

                timeStep++;
        }
}

// ----------------- Private methods -------------------

void RANDOM_BOOLEAN_NETWORK::Compute_Modularity(void) {

	double numArmArmEdges   = 0.0;
        double numHandHandEdges = 0.0;
        double numArmHandEdges  = 0.0;
        double numHandArmEdges  = 0.0;

        double numPossibleArmArmEdges   = (3.0 * NUM_NODES_PER_BODY_PART) * (3.0 * NUM_NODES_PER_BODY_PART);

        double numPossibleHandHandEdges = (4.0 * NUM_NODES_PER_BODY_PART) * (4.0 * NUM_NODES_PER_BODY_PART); 

        double numPossibleArmHandEdges  = (3.0 * NUM_NODES_PER_BODY_PART) * (4.0 * NUM_NODES_PER_BODY_PART);
 
        double numPossibleHandArmEdges  = (4.0 * NUM_NODES_PER_BODY_PART) * (3.0 * NUM_NODES_PER_BODY_PART);

	for (int e=0 ; e<numEdges; e++) {

		int originatingNode = edgeList->Get(e,0);
                int targetNode      = edgeList->Get(e,1);

		int originatingNodeIsInArm = originatingNode < (3*NUM_NODES_PER_BODY_PART);

                int targetNodeIsInArm      = targetNode      < (3*NUM_NODES_PER_BODY_PART);

		        if ( (originatingNodeIsInArm==true)  && (targetNodeIsInArm==true) )

			numArmArmEdges++;

                else    if ( (originatingNodeIsInArm==false) && (targetNodeIsInArm==false) )

                        numHandHandEdges++;

                else    if ( (originatingNodeIsInArm==true)  && (targetNodeIsInArm==false) )

                        numArmHandEdges++;

                else // if ( (originatingNodeIsInArm==false) && (targetNodeIsInArm==true) )

                        numHandArmEdges++;
	}

        double densityArmArmEdges   = numArmArmEdges   / numPossibleArmArmEdges;
        double densityHandHandEdges = numHandHandEdges / numPossibleHandHandEdges; 
        double densityArmHandEdges  = numArmHandEdges  / numPossibleArmHandEdges;
        double densityHandArmEdges  = numHandArmEdges  / numPossibleHandArmEdges;

	double numerator   = densityArmArmEdges + densityHandHandEdges;

	double denominator = densityArmHandEdges + densityHandArmEdges; 

	modularity = (1.0 + numerator) / (1.0 + denominator);
}

void RANDOM_BOOLEAN_NETWORK::Connection_Add(void) {

	int i = rand() % numNodes;
	int j = rand() % numNodes;

	while ( EdgeList_In(j,i) ) {

        	i = rand() % numNodes;
        	j = rand() % numNodes;
	}

	if ( (rand()%2) == 0 )

		EdgeList_AddEdge(j,i,+1.0);
	else
		EdgeList_AddEdge(j,i,-1.0);
}

void RANDOM_BOOLEAN_NETWORK::Connection_Add(int u) {

	int j = rand() % numNodes;

        while ( EdgeList_In(j,u) ) {

                j = rand() % numNodes; 
        }

        if ( (rand()%2)==0 ) {

                EdgeList_AddEdge(j,u,+1.0);
        }
        else {
                EdgeList_AddEdge(j,u,-1.0);
        }
}

void RANDOM_BOOLEAN_NETWORK::Connection_Remove(int u) {

        int j = rand() % numNodes;

        while ( !EdgeList_In(j,u) ) {

                j = rand() % numNodes;
        }

        EdgeList_Remove(j,u);
}

void RANDOM_BOOLEAN_NETWORK::EdgeList_AddEdge(int j, int i, double val) {

	edgeList->Set(numEdges,0,j);
	edgeList->Set(numEdges,1,i);
	edgeList->Set(numEdges,2,val);

	numEdges++;
}

int  RANDOM_BOOLEAN_NETWORK::EdgeList_In(int i, int j) {

        int found = false;
        int currEdge = 0;

        while ( (!found) && (currEdge<numEdges) ) {

		if (	(edgeList->vals[currEdge*edgeList->width]==i) &&
			(edgeList->vals[currEdge*edgeList->width+1]==j) )
                        found = true;
                else
                        currEdge++;
        }

        return( found );
}

void RANDOM_BOOLEAN_NETWORK::EdgeList_Initialize(void) {

        numEdges = 0;

        edgeList = new MATRIX(numNodes*numNodes,3,0.0);

        for (int i=0;i<START_NUM_CONNECTIONS;i++)

		Connection_Add();

	Compute_Modularity();
}

void RANDOM_BOOLEAN_NETWORK::EdgeList_Print(void) {

        for (int i=0;i<numEdges;i++) {

		printf("%d    \t",i);
                printf("%1.0f \t",edgeList->Get(i,0));
                printf("%1.0f \t",edgeList->Get(i,1));
                printf("%1.0f \n",edgeList->Get(i,2));

        }
}

void RANDOM_BOOLEAN_NETWORK::EdgeList_RegulateGenes(int timeStep) {

        for (int j=0;j<NUM_NODES;j++)

		nodeValues->vals[timeStep*NUM_NODES + j] = 0.0;

        for (int i=0;i<numEdges;i++) {

                int source    = int(edgeList->vals[i*edgeList->width + 0]);
                int target    = int(edgeList->vals[i*edgeList->width + 1]);

                double influence =  edgeList->vals[i*edgeList->width + 2];

                double nodeValue = nodeValues->vals[(timeStep-1)*NUM_NODES + source];

                double amountOfRegulation = nodeValue * influence;

                nodeValues->vals[timeStep*NUM_NODES + target] =
                nodeValues->vals[timeStep*NUM_NODES + target] +
                amountOfRegulation;
        }

        for (int j=0;j<NUM_NODES;j++)

		if ( nodeValues->vals[timeStep*NUM_NODES + j] > 0 )

			nodeValues->vals[timeStep*NUM_NODES + j] = +1.0;
		else
			nodeValues->vals[timeStep*NUM_NODES + j] = -1.0;
}

void RANDOM_BOOLEAN_NETWORK::EdgeList_Remove(int i, int j) {

        int found = false;
        int edgeIndex = 0;

        while ( !found ) {
                if (    (edgeList->Get(edgeIndex,0)==i) &&
                        (edgeList->Get(edgeIndex,1)==j) )
                        found = true;
                else
                        edgeIndex++;
        }

        for (int i=edgeIndex;i<(numEdges-1);i++) {

                edgeList->CopyRow(i+1,i);
        }
        numEdges--;
}

int RANDOM_BOOLEAN_NETWORK::Get_Angle_Index_2(int bodyPartIndex) {

        // Gray coding is employed

        int startNode = bodyPartIndex * NUM_NODES_PER_BODY_PART;

        if (            (nodeValues->Get(NUM_UPDATES-1,startNode)   <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1) <0) )
                return 0;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)   <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1)==1) )
                return 1;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)  ==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1)==1) )
                return 2;

        else // if (    (nodeValues->Get(NUM_UPDATES-1,startNode)  ==1) &&
             //         (nodeValues->Get(NUM_UPDATES-1,startNode+1) <0) )
                return 3;

}

int RANDOM_BOOLEAN_NETWORK::Get_Angle_Index_4(int bodyPartIndex) {

        // Gray coding is employed

        int startNode = bodyPartIndex * NUM_NODES_PER_BODY_PART;

        if (            (nodeValues->Get(NUM_UPDATES-1,startNode)   <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3) <0) )
                return 0;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)   <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3)==1) )
                return 1;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)   <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3)==1) )
                return 2;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)   <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3) <0) )
                return 3;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)   <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3) <0) )
                return 4;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)   <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3)==1) )
                return 5;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)   <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3)==1) )
                return 6;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)   <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3) <0) )
                return 7;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)  ==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3) <0) )
                return 8;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)  ==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3)==1) )
                return 9;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)  ==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3)==1) )
                return 10;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)  ==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3) <0) )
                return 11;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)  ==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3) <0) )
                return 12;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)  ==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2)==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3)==1) )
                return 13;

        else if (       (nodeValues->Get(NUM_UPDATES-1,startNode)  ==1) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+1) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+2) <0) &&
                        (nodeValues->Get(NUM_UPDATES-1,startNode+3)==1) )
                return 14;

        else // 1 0 0 0
                return 15;
}

void RANDOM_BOOLEAN_NETWORK::MutateGene(int u) {

        double N = numNodes;

        double ru = 0.0;

        for (int j=0;j<numNodes;j++)
                if ( EdgeList_In(j,u) )
                        ru++;

        double pu = (4*ru) / (4*ru + (N-ru) ); // Specialization paper, Eq. (5)

	if ( (float(rand())/float(RAND_MAX)) < pu )

                Connection_Remove(u);

        if ( (float(rand())/float(RAND_MAX)) < (1.0-pu) )

                Connection_Add(u);
}

void RANDOM_BOOLEAN_NETWORK::NodeValues_Initialize(void) {

	nodeValues = new MATRIX(NUM_UPDATES , NUM_NODES , 0.0);
}

#endif
