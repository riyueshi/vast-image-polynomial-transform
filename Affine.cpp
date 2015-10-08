//////////////////////////////////////////////////////////////////////
// Polynomial.cpp: implementation of the CPolynomial class.
// Author: Miller Zhang
// Email : imriyueshi@163.com
//////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <iostream>

#include "Affine.h"
#include "math.h"
#include "stdafx.h"
#include "opencv/cv.h"
#include "opencv2/highgui/highgui.hpp"

using namespace std;
using cv::Mat;

CAffine::CAffine()
{
}

CAffine::~CAffine()
{
}
bool CAffine::Compmuterpara(vector<Point2d_csu> &m_Arraylimage,vector<Point2d_csu> &m_Arrayrimage)
{
	if (m_Arraylimage.size() != m_Arrayrimage.size())
	{
		std::cout<<"error : size of points doesn't equal"<<std::endl;
		return false;
	}
	int m_nNumofPoint = m_Arraylimage.size();
	if (m_nNumofPoint < 3)
	{
		std::cout<<"too few points can be used, at least 3 point is requared"<<std::endl;
		return false;
	}
	Mat A = Mat::zeros(2*m_nNumofPoint,6,CV_32F);
	Mat l = Mat::zeros(2*m_nNumofPoint,1,CV_32F);
	for (int i = 0; i < m_nNumofPoint; i++)
	{
		A.at<float>(2*i,0) = 1;
		A.at<float>(2*i,1) = m_Arraylimage.at(i).x;
		A.at<float>(2*i,2) = m_Arraylimage.at(i).y;

		l.at<float>(2*i,0) = m_Arrayrimage.at(i).x;

		A.at<float>(2*i+1,3) = 1;
		A.at<float>(2*i+1,4) = m_Arraylimage.at(i).x;
		A.at<float>(2*i+1,5) = m_Arraylimage.at(i).y;

		l.at<float>(2*i+1,0) = m_Arrayrimage.at(i).y;

	}
	//cout<<"A: "<<A<<endl;
	//cout<<"l: "<<l<<endl;

	Mat params = (A.t()*A).inv()*A.t()*l;

	a0 = params.at<float>(0,0);
	a1 = params.at<float>(1,0);
	a2 = params.at<float>(2,0);

	b0 = params.at<float>(3,0);
	b1 = params.at<float>(4,0);
	b2 = params.at<float>(5,0);

	//cout<<"params: "<<params<<endl;
	return true;
}

//compute the source image pos (x,y) using the params which obtained by "Compmuterpara"	function
Point2d_csu CAffine::positiveMS(Point2d_csu ptsorce)
{
	Point2d_csu pt;	
	int xs = ptsorce.x;
	int ys = ptsorce.y;
	pt.x = a0 + a1*xs+a2*ys;
	pt.y = b0 + b1*xs+b2*ys;
	return pt;
}