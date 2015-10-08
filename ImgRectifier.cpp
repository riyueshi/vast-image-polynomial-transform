// ImgRectifier.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <string>
#include <vector>
#include <Windows.h>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <iomanip>

#include "Polynomial.h"
#include "opencv/cv.h"
#include "opencv2/highgui/highgui.hpp"

#include "gdal_priv.h"
#include "gdal_alg.h"

using namespace std;
using cv::Mat;


bool ImgPolyRectify(const string strRef, const string strReg, vector<Point2d_csu >RefPoint, vector<Point2d_csu>RegPoint, string out_image_name)
{
	GDALDataset *src_image;
	GDALDataset *ref_image;

	GDALAllRegister();

	src_image = (GDALDataset *) GDALOpen(strReg.c_str(), GA_ReadOnly );
	ref_image = (GDALDataset *) GDALOpen(strRef.c_str(), GA_ReadOnly );

	double transfer_params[6] = {0};
	ref_image->GetGeoTransform(transfer_params);

	CPolynomial param;
	CPolynomial inv_param;
	double g_min_x;
	double g_min_y;
	double g_max_x;
	double g_max_y;

	//compute the inv_param to get the range of result image
	if (!inv_param.Compmuterpara(RegPoint,RefPoint))
	{
		cout<<"inv_param compute failed!"<<endl;
		getchar();
		return false;
	}
	else
	{
		vector<Point2d_csu> src_point;
		src_point.resize(4);
		src_point.at(0) = Point2d_csu(0,0);
		src_point.at(1) = Point2d_csu(src_image->GetRasterXSize(),src_image->GetRasterYSize());
		src_point.at(2) = Point2d_csu(0,src_image->GetRasterYSize());
		src_point.at(3) = Point2d_csu(src_image->GetRasterXSize(),0);

		vector<Point2d_csu> tar_point;
		tar_point.resize(4);
		for (int i = 0; i < 4; i++)
		{					 
			tar_point.at(i) = inv_param.positiveMS(src_point.at(i));
		}
		vector<double> x_vec;
		vector<double> y_vec;
		for (int i = 0; i < 4; i++)
		{
			x_vec.push_back(tar_point.at(i).x);
			y_vec.push_back(tar_point.at(i).y);
		}
		g_min_x = *min_element(x_vec.begin(),x_vec.end());
		g_min_y = *min_element(y_vec.begin(),y_vec.end());
		g_max_x = *max_element(x_vec.begin(),x_vec.end());
		g_max_y = *max_element(y_vec.begin(),y_vec.end());
	}


	//compute the transfer params
	double x_value = transfer_params[0] + int(g_min_x)*transfer_params[1];
	double y_value = transfer_params[3] + int(g_min_y)*transfer_params[5];

	transfer_params[0] = x_value;
	transfer_params[3] = y_value;


	if (!param.Compmuterpara(RefPoint,RegPoint))
	{
		cout<<"params compute failed!"<<endl;
		getchar();
		return false;
	}
	else
	{
		vector<double> res_x;
		vector<double> res_y;
		vector<Point2d_csu> transed_points;
		for (int i = 0; i < RefPoint.size(); i++)
		{
			Point2d_csu src_pos = param.positiveMS(RefPoint.at(i));
			transed_points.push_back(src_pos);
			res_x.push_back(src_pos.x-RegPoint.at(i).x);
			res_y.push_back(src_pos.y-RegPoint.at(i).y);
		}
		double sum_x = std::accumulate(res_x.begin(), res_x.end(), 0.0);
		double mean_x = sum_x / res_x.size();
		double sq_sum_x = std::inner_product(res_x.begin(), res_x.end(), res_x.begin(), 0.0);
		double stdev_x = std::sqrt(sq_sum_x / res_x.size() - mean_x * mean_x);

		double sum_y = std::accumulate(res_y.begin(), res_y.end(), 0.0);
		double mean_y = sum_y / res_y.size();
		double sq_sum_y = std::inner_product(res_y.begin(), res_y.end(), res_y.begin(), 0.0);
		double stdev_y = std::sqrt(sq_sum_y / res_y.size() - mean_y * mean_y);	


		//output the transform infomation 
		ofstream res_writer("sigma.txt");
		res_writer<<setw(18)<<"sigma x: "<<stdev_x<<"   sigma y:"<<stdev_y<<endl;
		res_writer<<setiosflags(ios::fixed)<<setprecision(3);
		res_writer<<setw(18)<<"Point ID"
			<<setw(18)<<"res_x"
			<<setw(18)<<"res_y"
			<<setw(18)<<"ori_point_x"
			<<setw(18)<<"ori_point_y"
			<<setw(18)<<"transed_point_x"
			<<setw(18)<<"transed_point_y"<<endl;
		for (int i = 0; i < RefPoint.size(); i++)
		{
			res_writer<<setw(18)<<setprecision(0)<<i
				<<setprecision(3)
				<<setw(18)<<res_x.at(i)
				<<setw(18)<<res_y.at(i)
				<<setw(18)<<RegPoint.at(i).x
				<<setw(18)<<RegPoint.at(i).y
				<<setw(18)<<transed_points.at(i).x
				<<setw(18)<<transed_points.at(i).y<<endl;
		}
		res_writer.close();
	}

	//the image size to segment the source image , to set the same size of source image
	//means process the image without segment which can be applied to small size image
	const int nx_bolck_size = 3000;
	const int ny_bolck_size = 3000;

	const int nx_block = (g_max_x-g_min_x)/nx_bolck_size + 1;
	const int ny_block = (g_max_y-g_min_y)/ny_bolck_size + 1;

	const int ref_img_width = g_max_x-g_min_x;
	const int ref_img_height = g_max_y-g_min_y;

	GDALDataset *tar_image;
	GDALDriver *poDriver;  
	string fomat="GTiff";
	poDriver = GetGDALDriverManager()->GetDriverByName(fomat.c_str());
	tar_image = poDriver->Create(out_image_name.c_str() ,ref_img_width,ref_img_height,3 ,GDT_Byte,0);
	for (int block_index_x = 0; block_index_x < nx_block; block_index_x++)
	{
		for (int block_index_y = 0; block_index_y < ny_block; block_index_y++)
		{

			const int nx_size = (block_index_x+1)*nx_bolck_size<ref_img_width?nx_bolck_size:(ref_img_width-block_index_x*nx_bolck_size);
			const int ny_size = (block_index_y+1)*ny_bolck_size<ref_img_height?ny_bolck_size:(ref_img_height-block_index_y*ny_bolck_size);

			vector<double> x;
			vector<double> y;
			Point2d_csu src_pos1 = param.positiveMS(Point2d_csu(block_index_x*nx_bolck_size+g_min_x,block_index_y*ny_bolck_size+g_min_y));
			Point2d_csu src_pos2 = param.positiveMS(Point2d_csu(block_index_x*nx_bolck_size+nx_size+g_min_x,block_index_y*ny_bolck_size+g_min_y));
			Point2d_csu src_pos3 = param.positiveMS(Point2d_csu(block_index_x*nx_bolck_size+g_min_x,block_index_y*ny_bolck_size+ny_size+g_min_y));
			Point2d_csu src_pos4 = param.positiveMS(Point2d_csu(block_index_x*nx_bolck_size+nx_size+g_min_x,block_index_y*ny_bolck_size+ny_size+g_min_y));

			x.push_back(src_pos1.x);
			x.push_back(src_pos2.x);
			x.push_back(src_pos3.x);
			x.push_back(src_pos4.x);

			y.push_back(src_pos1.y);
			y.push_back(src_pos2.y);
			y.push_back(src_pos3.y);
			y.push_back(src_pos4.y);

			const int buffer_size(10);
			int minx = *min_element(x.begin(),x.end())-buffer_size;
			int miny = *min_element(y.begin(),y.end())-buffer_size;
			int maxx = *max_element(x.begin(),x.end())+buffer_size;
			int maxy = *max_element(y.begin(),y.end())+buffer_size;

			minx = minx>0?minx:0;
			miny = miny>0?miny:0;
			maxx = maxx<src_image->GetRasterXSize()?maxx:src_image->GetRasterXSize();
			maxy = maxy<src_image->GetRasterYSize()?maxy:src_image->GetRasterYSize();

			BYTE *src_image_data = new BYTE[(maxx-minx)*(maxy-miny)*3];
			BYTE *tar_image_data = new BYTE[nx_size*ny_size*3];

			int panBandMap [3]={3,2,1};
			src_image->RasterIO(GF_Read,minx,miny,maxx-minx,maxy-miny,src_image_data,maxx-minx,maxy-miny,GDT_Byte,3,panBandMap,3,(maxx-minx)*3,1);

			Mat temp_src = Mat((maxy-miny), (maxx-minx), CV_8UC3, src_image_data);
			Mat temp_tar = Mat::zeros(ny_size, nx_size, CV_8UC3);

			Mat mapx = Mat::zeros(temp_tar.rows,temp_tar.cols,CV_32FC1);
			Mat mapy = Mat::zeros(temp_tar.rows,temp_tar.cols,CV_32FC1);

			for (int row_index = 0; row_index < temp_tar.rows; row_index++)
			{
				for (int col_index = 0; col_index < temp_tar.cols; col_index++)
				{
					Point2d_csu src_pos = param.positiveMS(Point2d_csu(block_index_x*nx_bolck_size+col_index+g_min_x,block_index_y*ny_bolck_size+row_index+g_min_y));
					if (src_pos.x>-1&&src_pos.x<(src_image->GetRasterXSize())&&src_pos.y>-1&&src_pos.y<(src_image->GetRasterYSize()))
					{
						mapx.at<float>(row_index,col_index) = src_pos.x-minx;
						mapy.at<float>(row_index,col_index) = src_pos.y-miny;
					}
				}
			}
			cv::remap(temp_src,temp_tar,mapx,mapy,CV_INTER_LINEAR);
			tar_image->RasterIO(GF_Write,block_index_x*nx_bolck_size,block_index_y*ny_bolck_size,nx_size,ny_size,temp_tar.data,nx_size,ny_size,GDT_Byte,3,panBandMap,3,nx_size*3,1);
			delete []src_image_data;
			delete []tar_image_data;
		}
	}
	tar_image->SetGeoTransform(transfer_params);
	GDALClose( (GDALDatasetH) tar_image );
	GDALClose( (GDALDatasetH) src_image );
	GDALClose( (GDALDatasetH) ref_image );

	return true;


}


//update the mean square error of a point vector
bool PolyResCal(const vector<Point2d_csu>RefPoint, const vector<Point2d_csu>RegPoint, vector<double>&Res)
{
	Res.clear();
	CPolynomial param;
	if (!param.Compmuterpara(RefPoint,RegPoint))
	{
		cout<<"params compute failed!"<<endl;
		getchar();
		return false;
	}
	else
	{
		vector<double> res_x;
		vector<double> res_y;
		vector<Point2d_csu> transed_points;
		for (int i = 0; i < RefPoint.size(); i++)
		{
			Point2d_csu src_pos = param.positiveMS(RefPoint.at(i));
			transed_points.push_back(src_pos);
			res_x.push_back(src_pos.x-RegPoint.at(i).x);
			res_y.push_back(src_pos.y-RegPoint.at(i).y);
		}
		Res.clear();
		for (int i = 0; i < res_x.size(); i++)
		{
			Res.push_back(sqrt(res_x.at(i)*res_x.at(i)+res_y.at(i)*res_y.at(i)));
		}

	}
	return true;
}


//int _tmain(int argc, _TCHAR* argv[])
//{
//	GDALDataset *src_image;
//	GDALDataset *ref_image;
//
//	GDALAllRegister();
//
//	src_image = (GDALDataset *) GDALOpen("E:\\projectsE\\UAV_DATA\\OrthoPhoto\\block_0.tif", GA_ReadOnly );
//	ref_image = (GDALDataset *) GDALOpen("F:\\projectsF\\ImgRectifier\\19\\ref.tif", GA_ReadOnly );
//
//	vector<Point2d_csu> src_point;
//	vector<Point2d_csu> tar_point;
//
//	src_point.resize(7);
//	tar_point.resize(7);
//
//	//src_point.at(0).x = 948;
//	//src_point.at(0).y = 2585;
//	//src_point.at(1).x = 3127;
//	//src_point.at(1).y = 1006;
//	//src_point.at(2).x = 577;
//	//src_point.at(2).y = 1142;
//	//src_point.at(3).x = 3011;
//	//src_point.at(3).y = 1975;
//	//src_point.at(4).x = 2130;
//	//src_point.at(4).y = 1565;
//	//src_point.at(5).x = 422;
//	//src_point.at(5).y = 282;
//	//src_point.at(6).x = 2001;
//	//src_point.at(6).y = 623;
//
//	//tar_point.at(0).x = 1026;
//	//tar_point.at(0).y = 807;
//	//tar_point.at(1).x = 2373;
//	//tar_point.at(1).y = 2331;
//	//tar_point.at(2).x = 2232;
//	//tar_point.at(2).y = 475;
//	//tar_point.at(3).x = 1629;
//	//tar_point.at(3).y = 2297;
//	//tar_point.at(4).x = 1922;
//	//tar_point.at(4).y = 1669;
//	//tar_point.at(5).x = 2943;
//	//tar_point.at(5).y = 331;
//	//tar_point.at(6).x = 2661;
//	//tar_point.at(6).y = 1542;
//
//	//WenChuan image data
//	src_point.at(0).x = 2972;
//	src_point.at(0).y = 8896;
//	src_point.at(1).x = 4531;
//	src_point.at(1).y = 5921;
//	src_point.at(2).x = 5059;
//	src_point.at(2).y = 3032;
//	src_point.at(3).x = 7852;
//	src_point.at(3).y = 2526;
//	src_point.at(4).x = 2808;
//	src_point.at(4).y = 5809;
//	src_point.at(5).x = 2891;
//	src_point.at(5).y = 6986;
//	src_point.at(6).x = 6076;
//	src_point.at(6).y = 4249;
//
//	tar_point.at(0).x = 3584;
//	tar_point.at(0).y = 11244;
//	tar_point.at(1).x = 5537;
//	tar_point.at(1).y = 7456;
//	tar_point.at(2).x = 6269;
//	tar_point.at(2).y = 3839;
//	tar_point.at(3).x = 9814;
//	tar_point.at(3).y = 3204;
//	tar_point.at(4).x = 3400;
//	tar_point.at(4).y = 7352;
//	tar_point.at(5).x = 3487;
//	tar_point.at(5).y = 8827;
//	tar_point.at(6).x = 7530;
//	tar_point.at(6).y = 5366;
//
//
//	CPolynomial param;
//	if (!param.Compmuterpara(src_point,tar_point))
//	{
//		cout<<"params compute failed!"<<endl;
//		getchar();
//		return -1;
//	}
//	const int nx_bolck_size = 2000;
//	const int ny_bolck_size = 2000;
//
//	const int nx_block = ref_image->GetRasterXSize()/nx_bolck_size + 1;
//	const int ny_block = ref_image->GetRasterYSize()/ny_bolck_size + 1;
//
//	const int ref_img_width = ref_image->GetRasterXSize();
//	const int ref_img_height = ref_image->GetRasterYSize();
//
//	GDALDataset *tar_image;
//	GDALDriver *poDriver;  
//	string fomat="GTiff";
//	poDriver = GetGDALDriverManager()->GetDriverByName(fomat.c_str());
//	tar_image = poDriver->Create("result.tif",ref_img_width,ref_img_height,3 ,GDT_Byte,0);
//	for (int block_index_x = 0; block_index_x < nx_block; block_index_x++)
//	{
//		for (int block_index_y = 0; block_index_y < ny_block; block_index_y++)
//		{
//
//			const int nx_size = (block_index_x+1)*nx_bolck_size<ref_img_width?nx_bolck_size:(ref_img_width-block_index_x*nx_bolck_size);
//			const int ny_size = (block_index_y+1)*ny_bolck_size<ref_img_height?ny_bolck_size:(ref_img_height-block_index_y*ny_bolck_size);
//
//			vector<double> x;
//			vector<double> y;
//			Point2d_csu src_pos1 = param.positiveMS(Point2d_csu(block_index_x*nx_bolck_size,block_index_y*ny_bolck_size));
//			Point2d_csu src_pos2 = param.positiveMS(Point2d_csu(block_index_x*nx_bolck_size+nx_size,block_index_y*ny_bolck_size));
//			Point2d_csu src_pos3 = param.positiveMS(Point2d_csu(block_index_x*nx_bolck_size,block_index_y*ny_bolck_size+ny_size));
//			Point2d_csu src_pos4 = param.positiveMS(Point2d_csu(block_index_x*nx_bolck_size+nx_size,block_index_y*ny_bolck_size+ny_size));
//
//			x.push_back(src_pos1.x);
//			x.push_back(src_pos2.x);
//			x.push_back(src_pos3.x);
//			x.push_back(src_pos4.x);
//
//			y.push_back(src_pos1.y);
//			y.push_back(src_pos2.y);
//			y.push_back(src_pos3.y);
//			y.push_back(src_pos4.y);
//
//			const int buffer_size(10);
//			int minx = *min_element(x.begin(),x.end())-buffer_size;
//			int miny = *min_element(y.begin(),y.end())-buffer_size;
//			int maxx = *max_element(x.begin(),x.end())+buffer_size;
//			int maxy = *max_element(y.begin(),y.end())+buffer_size;
//
//			minx = minx>0?minx:0;
//			miny = miny>0?miny:0;
//			maxx = maxx<src_image->GetRasterXSize()?maxx:src_image->GetRasterXSize();
//			maxy = maxy<src_image->GetRasterYSize()?maxy:src_image->GetRasterYSize();
//
//			BYTE *src_image_data = new BYTE[(maxx-minx)*(maxy-miny)*3];
//			BYTE *tar_image_data = new BYTE[nx_size*ny_size*3];
//
//			int panBandMap [3]={3,2,1};
//			src_image->RasterIO(GF_Read,minx,miny,maxx-minx,maxy-miny,src_image_data,maxx-minx,maxy-miny,GDT_Byte,3,panBandMap,3,(maxx-minx)*3,1);
//
//			Mat temp_src = Mat((maxy-miny), (maxx-minx), CV_8UC3, src_image_data);
//			Mat temp_tar = Mat(ny_size, nx_size, CV_8UC3);
//
//			Mat mapx(temp_tar.rows,temp_tar.cols,CV_32FC1);
//			Mat mapy(temp_tar.rows,temp_tar.cols,CV_32FC1);
//
//			for (int row_index = 0; row_index < temp_tar.rows; row_index++)
//			{
//				for (int col_index = 0; col_index < temp_tar.cols; col_index++)
//				{
//					Point2d_csu src_pos = param.positiveMS(Point2d_csu(block_index_x*nx_bolck_size+col_index,block_index_y*ny_bolck_size+row_index));
//					if (src_pos.x>0&&src_pos.x<(src_image->GetRasterXSize())&&src_pos.y>0&&src_pos.y<(src_image->GetRasterYSize()))
//					{
//						mapx.at<float>(row_index,col_index) = src_pos.x-minx;
//						mapy.at<float>(row_index,col_index) = src_pos.y-miny;
//					}
//				}
//			}
//			cv::remap(temp_src,temp_tar,mapx,mapy,CV_INTER_LINEAR);
//			tar_image->RasterIO(GF_Write,block_index_x*nx_bolck_size,block_index_y*ny_bolck_size,nx_size,ny_size,temp_tar.data,nx_size,ny_size,GDT_Byte,3,panBandMap,3,nx_size*3,1);
//			delete []src_image_data;
//			delete []tar_image_data;
//		}
//	}
//	GDALClose( (GDALDatasetH) tar_image );
//	GDALClose( (GDALDatasetH) src_image );
//	GDALClose( (GDALDatasetH) ref_image );
//
//	return 0;
//}

