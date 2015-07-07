#include "matrix.h"
#include "object.h"

#ifndef _ROBOT_H
#define _ROBOT_H

class ROBOT {

public:
	int numPossibleAngles;
	MATRIX *robot;

public:
	ROBOT(void);
	~ROBOT(void);
	int    Claws_Cross(void);
	void   ComputePositions(void);
	double Get_Diff_Between_Fingertip_Dist_And_Circumference_Of(OBJECT *object);
	double Get_LeftFingerTip_Position_X(void);
        double Get_LeftFingerTip_Position_Y(void);
	double Get_LeftFingerTip_Distance_From_Circumference_Of(OBJECT *object);
        double Get_RightFingerTip_Position_X(void);
        double Get_RightFingerTip_Position_Y(void);
        double Get_RightFingerTip_Distance_From_Circumference_Of(OBJECT *object);
	void   Print(void);
	void   Save(char *fileName);
	void   SetAngle(int bodyPartIndex, double angle);
	void   SetLength(int bodyPartIndex, double length);
	void   UpdateAngles(void);

private:
	int    Line_Segments_Intersect(int segment1, int segment2);
};

#endif
