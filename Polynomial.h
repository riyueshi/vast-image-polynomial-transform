// Polynomial.h: interface for the CPolynomial class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include <string>
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

class CPolynomial  
{
public:
	CPolynomial();
	virtual ~CPolynomial();
	bool Compmuterpara(const std::vector<Point2d_csu> &m_Arraylimage,const std::vector<Point2d_csu> &m_Arrayrimage);
	Point2d_csu CPolynomial::positiveMS(Point2d_csu ptsorce);
public:
	double a0,a1,a2,a3,a4,a5;
	double b0,b1,b2,b3,b4,b5;

};

#ifndef POLYRECTIFY_DLL_H
#define POLYRECTIFY_DLL_H

extern "C" __declspec(dllexport) bool ImgPolyRectify(const std::string strRef, const std::string strReg, const std::vector<Point2d_csu >RefPoint, const std::vector<Point2d_csu>RegPoint, std::string out_image_name);

extern "C" __declspec(dllexport) bool PolyResCal(const std::vector<Point2d_csu>RefPoint, const std::vector<Point2d_csu>RegPoint, std::vector<double>&Res);

#endif