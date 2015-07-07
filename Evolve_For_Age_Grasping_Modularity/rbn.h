#include "matrix.h"
#include "robot.h"

#ifndef _RANDOM_BOOLEAN_NETWORK_H
#define _RANDOM_BOOLEAN_NETWORK_H

class RANDOM_BOOLEAN_NETWORK {

public:
        MATRIX *nodeValues;

private:
	int numNodes;
        int numEdges;
        MATRIX *edgeList;
	double modularity;

public:
	RANDOM_BOOLEAN_NETWORK(void);
	~RANDOM_BOOLEAN_NETWORK(void);
	void CopyFrom( RANDOM_BOOLEAN_NETWORK *parent );
	double Distance_From_Target_Pattern(int objectIndex);
	double Fraction_Of_Nodes_Unchanging(void);
	int  Get_Angle_Index(int bodyPartIndex);
	int  Get_NodeValue(int nodeIndex);
	int  In_Point_Attractor(void);
	double Get_Modularity(void);
        void Mutate(void);
	void Perturb_Node(int nodeIndex);
	void Print(void);
	void Set_Initial_Conditions(int environmentIndex);
	void Update(void);

private:
	void Compute_Modularity(void);
	void Connection_Add(void);
        void Connection_Add(int u);
        void Connection_Remove(int u);
        void EdgeList_AddEdge(int j, int i, double val);
        int  EdgeList_In(int i, int j);
        void EdgeList_Initialize(void);
        void EdgeList_Print(void);
        void EdgeList_RegulateGenes(int timeStep);
        void EdgeList_Remove(int i, int j);
	int  Get_Angle_Index_2(int bodyPartIndex);
	int  Get_Angle_Index_4(int bodyPartIndex);
        void MutateGene(int u);
	void NodeValues_Initialize(void);
};

#endif
