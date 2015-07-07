#ifndef _OBJECT_H
#define _OBJECT_H

class OBJECT {

private:
	double x;
	double y;
	double radius;

public:
	OBJECT(int objectIndex);
	~OBJECT(void);
	double Get_X(void);
	double Get_Y(void);
	double Get_Radius(void);
};

#endif
