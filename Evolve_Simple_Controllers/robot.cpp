/* ---------------------------------------------------
   FILE:     bodyPlan.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#ifndef _ROBOT_CPP
#define _ROBOT_CPP

#include "math.h"

#include "constants.h"
#include "robot.h"
#include "matrix.h"

ROBOT::ROBOT(void) {

	robot = new MATRIX(NUM_BODY_PARTS,5,0.0); 

	// Parent object    length               angle               x                  y
	robot->Set(0,0,-1); robot->Set(0,1,1);   robot->Set(0,2,0);  robot->Set(0,3,0); robot->Set(0,4,0);
        robot->Set(1,0,0);  robot->Set(1,1,1);   robot->Set(1,2,0);  robot->Set(1,3,0); robot->Set(1,4,0);
        robot->Set(2,0,1);  robot->Set(2,1,1);   robot->Set(2,2,0);  robot->Set(2,3,0); robot->Set(2,4,0);
	// Arm

        robot->Set(3,0,2);  robot->Set(3,1,0.5); robot->Set(3,2,0);  robot->Set(3,3,0); robot->Set(3,4,0); // left
        robot->Set(4,0,2);  robot->Set(4,1,0.5); robot->Set(4,2,0);  robot->Set(4,3,0); robot->Set(4,4,0); // right
	// Proximate phalanges

        robot->Set(5,0,3);  robot->Set(5,1,0.5); robot->Set(5,2,0);  robot->Set(5,3,0); robot->Set(5,4,0); // left
        robot->Set(6,0,4);  robot->Set(6,1,0.5); robot->Set(6,2,0);  robot->Set(6,3,0); robot->Set(6,4,0); // right
	// Distal phalanges
}

ROBOT::~ROBOT(void) {

	if ( robot ) {
		delete robot;
		robot = NULL;
	}
}

int  ROBOT::Claws_Cross(void) {

	int ld_rp = Line_Segments_Intersect(5,4); // Left distal and right proximal phalanges intersect.
        int ld_rd = Line_Segments_Intersect(5,6); // Left distal and right distal phalanges intersect.


        int rd_lp = Line_Segments_Intersect(6,4); // Right distal and left proximal phalanges intersect.

	return ( ld_rp || ld_rd || rd_lp );
}

void ROBOT::ComputePositions(void) {

	double length, angle;
	double x,y;
	length = robot->Get(0,1);
	angle  = robot->Get(0,2);

	x = length * sin(3.14159*angle/180.0);
	y = length * cos(3.14159*angle/180.0);

	robot->Set(0,3,x);
	robot->Set(0,4,y);

	for (int i=1;i<NUM_BODY_PARTS;i++) {
		int parent = int(robot->Get(i,0));
		x = robot->Get(parent,3);
		y = robot->Get(parent,4);
		length = robot->Get(i,1);
		angle  = robot->Get(i,2);
	
		x = x + length * sin(3.14159*angle/180.0);
		y = y + length * cos(3.14159*angle/180.0);
		robot->Set(i,3,x);
		robot->Set(i,4,y);
	}
}

double ROBOT::Get_Diff_Between_Fingertip_Dist_And_Circumference_Of(OBJECT *object) {

        double xLeft = Get_LeftFingerTip_Position_X();
        double yLeft = Get_LeftFingerTip_Position_Y();

        double xRight = Get_RightFingerTip_Position_X();
        double yRight = Get_RightFingerTip_Position_Y();

	double fingertipDist = sqrt( pow(xLeft-xRight,2.0) + pow(yLeft-yRight,2.0) );

	double objRadius = object->Get_Radius();

	double objCircumference = objRadius * 2.0;

	double diff = fabs( fingertipDist - objCircumference );

	return diff;
}

double ROBOT::Get_LeftFingerTip_Distance_From_Circumference_Of(OBJECT *object) {

	double x = Get_LeftFingerTip_Position_X();
	double y = Get_LeftFingerTip_Position_Y();

	double objX = object->Get_X();
	double objY = object->Get_Y();

	double distFromCenter = sqrt( pow(x-objX,2.0) + pow(y-objY,2.0) );
	
	double distFromCircumference = fabs( distFromCenter - object->Get_Radius() );

	return distFromCircumference;
}

double ROBOT::Get_LeftFingerTip_Position_X(void) {

	return( robot->Get(5,3) );
}

double ROBOT::Get_LeftFingerTip_Position_Y(void) {

        return( robot->Get(5,4) );
}

double ROBOT::Get_RightFingerTip_Distance_From_Circumference_Of(OBJECT *object) {

        double x = Get_RightFingerTip_Position_X();
        double y = Get_RightFingerTip_Position_Y();

        double objX = object->Get_X();
        double objY = object->Get_Y();

        double distFromCenter = sqrt( pow(x-objX,2.0) + pow(y-objY,2.0) );

        double distFromCircumference = fabs( distFromCenter - object->Get_Radius() );

        return distFromCircumference;
}

double ROBOT::Get_RightFingerTip_Position_X(void) {

        return( robot->Get(6,3) );
}

double ROBOT::Get_RightFingerTip_Position_Y(void) {

        return( robot->Get(6,4) );
}

void ROBOT::Print(void) {

	robot->Print(3);
}

void ROBOT::Save(char *fileName) {

	ofstream *outFile = new ofstream(fileName,ios::trunc);

	robot->WriteWithoutPreamble(outFile);

	(*outFile) << "\n";

	outFile->close();
	delete outFile;
	outFile = NULL;
}

void ROBOT::SetAngle(int bodyPartIndex, double angle) {

	robot->Set(bodyPartIndex,2,angle);
}

void ROBOT::SetLength(int bodyPartIndex, double length) {

	robot->Set(bodyPartIndex,1,length);
}

void ROBOT::UpdateAngles(void) {

	for (int i=1;i<NUM_BODY_PARTS;i++) {

		int parent 		= int(robot->Get(i,0));
		double parentsAngle 	= robot->Get(parent,2);

		robot->Add(i,2,parentsAngle);	
	}
}

// ------------------- Private methods -------------------------

int ROBOT::Line_Segments_Intersect(int segment1Index, int segment2Index) {

	int parent1Index = robot->Get(segment1Index,0);
        int parent2Index = robot->Get(segment2Index,0);

	double Ax = robot->Get(parent1Index,3);
        double Ay = robot->Get(parent1Index,4);

        double Bx = robot->Get(segment1Index,3);
        double By = robot->Get(segment1Index,4);

        double Cx = robot->Get(parent2Index,3);
        double Cy = robot->Get(parent2Index,4);

        double Dx = robot->Get(segment2Index,3);
        double Dy = robot->Get(segment2Index,4);

	double  distAB, theCos, theSin, newX, ABpos ;

	//  Fail if either line segment is zero-length.
	if ( ((Ax==Bx) && (Ay==By)) || ((Cx==Dx) && (Cy==Dy)) )

		return false; 

	//  Fail if the segments share an end-point.
	if ( ((Ax==Cx) && (Ay==Cy)) || ((Bx==Cx) && (By==Cy))
	 ||  ((Ax==Dx) && (Ay==Dy)) || ((Bx==Dx) && (By==Dy)) )

    		return false; 

	//  (1) Translate the system so that point A is on the origin.
	Bx-=Ax; By-=Ay;
	Cx-=Ax; Cy-=Ay;
	Dx-=Ax; Dy-=Ay;

	//  Discover the length of segment A-B.
	distAB=sqrt(Bx*Bx+By*By);

	//  (2) Rotate the system so that point B is on the positive X axis.
	theCos=Bx/distAB;
	theSin=By/distAB;
	newX=Cx*theCos+Cy*theSin;
	Cy  =Cy*theCos-Cx*theSin; Cx=newX;
	newX=Dx*theCos+Dy*theSin;
	Dy  =Dy*theCos-Dx*theSin; Dx=newX;

  	//  Fail if segment C-D doesn't cross line A-B.
  	if ( ((Cy<0.) && (Dy<0.)) || ((Cy>=0.) && (Dy>=0.)) )
		
		return false;

	//  (3) Discover the position of the intersection point along line A-B.
	ABpos=Dx+(Cx-Dx)*Dy/(Dy-Cy);

	//  Fail if segment C-D crosses line A-B outside of segment A-B.
	if ( (ABpos<0.) || (ABpos>distAB) )

		return false;

	// Success
	return true;
}


#endif
