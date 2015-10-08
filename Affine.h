#pragma once

#include <vector>

struct Point2d_csu
{
	Point2d_csu(){};
	Point2d_csu(double _x, double _y)
	{
		x = _x;
		y = _y;
	}
	double x;
	double y;
};

class CAffine  
{
public:
	CAffine();
	virtual ~CAffine();
	bool Compmuterpara(std::vector<Point2d_csu> &m_Arraylimage,std::vector<Point2d_csu> &m_Arrayrimage);
	Point2d_csu CAffine::positiveMS(Point2d_csu ptsorce);
public:
	double a0,a1,a2;
	double b0,b1,b2;

};

