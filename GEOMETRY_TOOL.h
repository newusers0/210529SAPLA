#pragma once
#ifndef GEOMETRY_TOOL_H
#define GEOMETRY_TOOL_H

#include "pch.h"

TEMPLATE
class GEOMETRY_TOOL
{

public:
	struct POINT;//1801178
	struct GEOMETRY::POINT center_geometry;//Function1
	struct GEOMETRY::POINT p0;//Function2

	struct ACCUMULATE_ID;
	struct ACCUMULATE_VALUE;
	//struct ID_INCREASE;
	struct ValueIncrease;//190607  min point value to max point value
	struct ValueDecrease;//190608  max point value to min point value

	template<typename T, typename S>
	S& getAverage(const T& const original_array, S& const result_argu);

	/*======================== Triangle ===================================*/
	//210402
	inline long double get_abs_area_triangle(const typename GEOMETRY::POINT& const point1, const  typename GEOMETRY::POINT& const point2, const  typename GEOMETRY::POINT& const point3);
	/*=====================================================================*/
	/*************************************************************/
	//Function1 Clockwise need center point
	static bool PointCmp(POINT a, POINT b);
	static bool PointCmpCenter(POINT& center, const POINT& a, const POINT& b);//190227
	static bool PointCmp(POINT& center, POINT a, POINT b);//190227
	struct CMP;//180123
	bool Cmp(const POINT& a, const POINT& b);
	POINT Gravity(const vector <GEOMETRY::POINT>& const p);// get center point
	//void sort(vector <POINT>& vPoints);
	void ClockwiseSortPoints(vector <GEOMETRY::POINT>& vPOINTs);
	////////////////////////////////////////////////////////////////

	/*************************************************************/
	//Function2 clock wise
	struct COMPARE_POINTS;
	static int get_quadrant(const POINT& p);
	static double get_clockwise_angle(const POINT& p);
	bool compare_points(const POINT& a, const POINT& b);
	//POINT center(0, 0);
	////////////////////////////////////////////////////////////////

	/*************************************************************/
	//Function3 Only compute area
	double ComputePolygonArea(const vector<POINT>& points);//190225
	///////////////////////////////////////////////////////////////

	/*******************************************************************/
	//Convex Hull (Graham Scan) 190226
	// A utility function to find next to top in a stack
	POINT nextToTop(stack<POINT>& S);
	// A utility function to swap two points 190226
	void swap(POINT& p1, POINT& p2);
	// A utility function to return square of distance, between p1 and p2 190226
	static double distSq(POINT p1, POINT p2);
	// To find orientation of ordered triplet (p, q, r). 190226
	// The function returns following values
	// 0 --> p, q and r are colinear
	// 1 --> Clockwise
	// 2 --> Counterclockwise
	static int orientation(POINT p, POINT q, POINT r);
	// A function used by library function qsort() to sort an array of 190226
	// points with respect to the first point
	int compare(const void* vp1, const void* vp2);
	// The same as up function compare();
	struct COMPARE_POLOR_POINTS;
	// Prints convex hull of a set of n points. 190227
	void convexHull(POINT points[], int n, std::vector <GEOMETRY::POINT>& points_convex);
	/////////////////////////////////////////////////////////////////////////

	/******************Convex or concave**********************************/
	//190305
	bool isHollow(std::vector<GEOMETRY::POINT>& curveloopPoints);
	//190306
	//int isConvex(POINT *p, int n);
	int isConvex(std::vector<GEOMETRY::POINT>& curveloopPoints);
	//////////////////////////////////////////////////////////////////////

	/*******************Area of Circle**************************************/
	//190312
	template<typename T>
	inline T area_of_a_circle(T r);
	//190308
	double get2PointsDist(const GEOMETRY::POINT& const point_1, const GEOMETRY::POINT& const point_2);
	//190308
	double getPointsMaxDist(const std::vector<GEOMETRY::POINT>& const points, const  GEOMETRY::POINT& const center_point);
	
	/////////////////////////////////////////////////////////////////////////

	/************************ 190311 Area of Circle******************************/
	
	int sgn(double x);
	double get_distance(const POINT a, const POINT b);
	POINT get_circle_center(const GEOMETRY::POINT a, const GEOMETRY::POINT b, const GEOMETRY::POINT c);
	
	void min_cover_circle(const std::vector<GEOMETRY::POINT>& const points, GEOMETRY::POINT& const center_point, double& const circle_radius);

};

#include "GEOMETRY_TOOL.cpp"

#endif

