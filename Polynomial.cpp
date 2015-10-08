//////////////////////////////////////////////////////////////////////
// Polynomial.cpp: implementation of the CPolynomial class.
// Author: Miller Zhang
// Email : imriyueshi@163.com
//////////////////////////////////////////////////////////////////////
#include <vector>
#include <string>
#include <iostream>

#include "stdafx.h"
#include "Polynomial.h"
#include "opencv/cv.h"
#include "opencv2/highgui/highgui.hpp"

using namespace std;
using cv::Mat;


CPolynomial::CPolynomial()
{

}

CPolynomial::~CPolynomial()
{

}

bool CPolynomial::Compmuterpara(const vector<Point2d_csu> &m_Arraylimage,const vector<Point2d_csu> &m_Arrayrimage)
{
	if (m_Arraylimage.size() != m_Arrayrimage.size())
	{
		std::cout<<"error : size of points doesn't equal"<<std::endl;
		return false;
	}
	int m_nNumofPoint = m_Arraylimage.size();
	if (m_nNumofPoint < 6)
	{
		std::cout<<"too few points can be used, at least 6 point is requared"<<std::endl;
		return false;
	}
	Mat A = Mat::zeros(2*m_nNumofPoint,12,CV_32F);
	Mat l = Mat::zeros(2*m_nNumofPoint,1,CV_32F);
	for (int i = 0; i < m_nNumofPoint; i++)
	{
		A.at<float>(2*i,0) = 1;
		A.at<float>(2*i,1) = m_Arraylimage.at(i).x * m_Arraylimage.at(i).x;
		A.at<float>(2*i,2) = m_Arraylimage.at(i).x;
		A.at<float>(2*i,3) = m_Arraylimage.at(i).x * m_Arraylimage.at(i).y;
		A.at<float>(2*i,4) = m_Arraylimage.at(i).y;
		A.at<float>(2*i,5) = m_Arraylimage.at(i).y * m_Arraylimage.at(i).y;

		l.at<float>(2*i,0) = m_Arrayrimage.at(i).x;

		A.at<float>(2*i+1,6) = 1;
		A.at<float>(2*i+1,7) = m_Arraylimage.at(i).x * m_Arraylimage.at(i).x;
		A.at<float>(2*i+1,8) = m_Arraylimage.at(i).x;
		A.at<float>(2*i+1,9) = m_Arraylimage.at(i).x * m_Arraylimage.at(i).y;
		A.at<float>(2*i+1,10) = m_Arraylimage.at(i).y;
		A.at<float>(2*i+1,11) = m_Arraylimage.at(i).y * m_Arraylimage.at(i).y;

		l.at<float>(2*i+1,0) = m_Arrayrimage.at(i).y;

	}

	Mat params = (A.t()*A).inv()*A.t()*l;

	a0 = params.at<float>(0,0);
	a1 = params.at<float>(1,0);
	a2 = params.at<float>(2,0);
	a3 = params.at<float>(3,0);
	a4 = params.at<float>(4,0);
	a5 = params.at<float>(5,0);

	b0 = params.at<float>(6,0);
	b1 = params.at<float>(7,0);
	b2 = params.at<float>(8,0);
	b3 = params.at<float>(9,0);
	b4 = params.at<float>(10,0);
	b5 = params.at<float>(11,0);

	return true;
}

//compute the source image pos (x,y) using the params which obtained by "Compmuterpara"	function
Point2d_csu CPolynomial::positiveMS(Point2d_csu ptsorce)
{
	Point2d_csu pt;	
	int xs = ptsorce.x;
	int ys = ptsorce.y;
	pt.x = a0 + a1*xs*xs+a2*xs+a3*xs*ys+a4*ys+a5*ys*ys;
	pt.y = b0 + b1*xs*xs+b2*xs+b3*xs*ys+b4*ys+b5*ys*ys;
	return pt;
}