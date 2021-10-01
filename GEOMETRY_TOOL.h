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
	//200206
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
	//190306
	//int isConvex(POINT *p, int n);
	//int isConvex(std::vector<GEOMETRY::POINT>& curveloopPoints);
	/////////////////////////////////////////////////////////////////////////

	/************************ 190311 Area of Circle******************************/
	//190311 https://www.cnblogs.com/Howe-Young/p/4477550.html
	int sgn(double x);
	double get_distance(const POINT a, const POINT b);//两点之间的距离
	POINT get_circle_center(const GEOMETRY::POINT a, const GEOMETRY::POINT b, const GEOMETRY::POINT c); //得到triangle外接圆的圆心
	//p表示定点, n表示顶点的个数, c代表最小覆盖圆圆心, r是半径
	void min_cover_circle(const std::vector<GEOMETRY::POINT>& const points, GEOMETRY::POINT& const center_point, double& const circle_radius);//找最小覆盖圆(这里没有用全局变量p[], 因为是为了封装一个函数便于调用)
	/////////////////////////////////////////////////////////////////////////////

};

//#include "GEOMETRY_TOOL.cpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**************************************************************************************************************************************************************/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//************************************
// Method:POINT
// Qualifier:
// Input:
// Output:
// date:190227
// author:
//************************************
TEMPLATE
struct GEOMETRY::POINT {//1801178
	double id = INF;
	double value = INF;
	//GEOMETRY::POINT* center;

	POINT(double x, double y) : id(x), value(y) {}
	POINT() : id(INF), value(INF) {}

	~POINT() {
		id = INF;
		value = INF;
	}
};

//************************************
// Method:ACCUMULATE_ID
// Qualifier:
// Input:
// Output:
// date:190307
// author:
//************************************
TEMPLATE
struct  GEOMETRY::ACCUMULATE_ID {
	double operator()(GEOMETRY::POINT& a, GEOMETRY::POINT& b) {
		return a.id + b.id;
	}
};

//************************************
// Method:ACCUMULATE_VALUE
// Qualifier:
// Input:
// Output:
// date:190307
// author:
//************************************
TEMPLATE
struct  GEOMETRY::ACCUMULATE_VALUE {
	double operator()(GEOMETRY::POINT a, GEOMETRY::POINT b) {
		return a.value + b.value;
	}
};

//************************************
// Method:ValueIncrease
// Qualifier:get min value of point
// Input:
// Output:
// date:190607
// author:
//************************************
TEMPLATE
struct  GEOMETRY::ValueIncrease {//190607  min point value to max point value
	bool operator()(GEOMETRY::POINT a, GEOMETRY::POINT b) {
		return a.value < b.value;
	}
};

//************************************
// Method:ValueDecrease
// Qualifier:get max value of point
// Input:
// Output:
// date:190608
// author:
//************************************
TEMPLATE
struct  GEOMETRY::ValueDecrease {//190608  max point value to min point value
	bool operator()(GEOMETRY::POINT a, GEOMETRY::POINT b)  const {
		return a.value > b.value;
	}
};

//************************************
// Method:PointCmpCenter
// Qualifier:
// Input:
// Output:
// date:190227
// author:
//************************************
TEMPLATE
template<typename T, typename S>
S& GEOMETRY::getAverage(const T& const original_array, S& const result_argu) {
	double id_sum = 0;
	double value_sum = 0;

	for (auto&& i : original_array) {
		id_sum += i.id;
		value_sum += i.value;
	}

	result_argu.id = id_sum / original_array.size();
	result_argu.value = value_sum / original_array.size();

	return result_argu;
}

//200206
//************************************
// Method:get_abs_area_triangle
// Qualifier://https://zhuanlan.zhihu.com/p/25793392
// Input:
// Output:
// date:200206
// author:
//************************************
TEMPLATE
inline long double GEOMETRY::get_abs_area_triangle(const typename GEOMETRY::POINT& const point1, const typename GEOMETRY::POINT& const point2, const  typename GEOMETRY::POINT& const point3) {
#ifdef _DEBUG
	assert(point1.id != INF && point3.id != INF && point2.value != INF && point3.value != INF && point1.value != INF && point3.value != INF && point2.id != INF && point3.id != INF);
#endif
	// 2 times of triangle abc area
	//double area_abc = (point1.id - point3.id) * (point2.value - point3.value) - (point1.value - point3.value) * (point2.id - point3.id);
	//return fabs(area_abc * 0.5);

	return fabs(((point1.id - point3.id) * (point2.value - point3.value) - (point1.value - point3.value) * (point2.id - point3.id)) * 0.5);
}

//************************************
// Method:PointCmpCenter
// Qualifier:
// Input:
// Output:
// date:190227
// author:
//************************************
TEMPLATE
bool GEOMETRY::PointCmp(POINT a, POINT b) {
	GEOMETRY::POINT center;
	//center.value = center1.value;
	//center.id = center1.id;
	//GEOMETRY::center1.value = 0;
	//GEOMETRY::center1.id = 0;
	if (a.id >= 0 && b.id < 0)
		return true;
	if (a.id == 0 && b.id == 0)
		return a.value > b.value;
	/////////////////////////// vector product of OA and OB
	int det = (a.id - center.id) * (b.value - center.value) - (b.id - center.id) * (a.value - center.value);
	if (det < 0)
		return true;
	if (det > 0)
		return false;
	/////////////////////////////////// OA and OB in the same line, according to the distance to decide priority.
	int d1 = (a.id - center.id) * (a.id - center.id) + (a.value - center.value) * (a.value - center.value);
	int d2 = (b.id - center.id) * (b.id - center.value) + (b.value - center.value) * (b.value - center.value);
	return d1 > d2;
}

//************************************
// Method:PointCmpCenter
// Qualifier:
// Input:
// Output:
// date:190227
// author:
//************************************
TEMPLATE
bool GEOMETRY::PointCmpCenter(POINT& center, const POINT& a, const POINT& b) {//190227
	//GEOMETRY::POINT center;
	assert(center.id != NULL && center.value != NULL);
	//assert(center.id == a.center->id && center.value == b.center->value);
	//assert(a.center->id == b.center->id && a.center->value == b.center->value);
	//center.id = a.center->id;
	//center.value = a.center->value;
	//cout << center.id << " " << center.value<<"::::::"<<a.id<<" "<<a.value<<" "<<b.id<<" "<<b.value << endl;
	//center.value = center1.value;
	//center.id = center1.id;
	//GEOMETRY::center1.value = 0;
	//GEOMETRY::center1.id = 0;
	if (a.id >= 0 && b.id < 0)
		return true;
	if (a.id == 0 && b.id == 0)
		return a.value > b.value;
	/////////////////////////// vector product of OA and OB
	int det = (a.id - center.id) * (b.value - center.value) - (b.id - center.id) * (a.value - center.value);
	if (det < 0)
		return true;
	if (det > 0)
		return false;
	/////////////////////////////////// OA and OB in the same line, according to the distance to decide priority.
	int d1 = (a.id - center.id) * (a.id - center.id) + (a.value - center.value) * (a.value - center.value);
	int d2 = (b.id - center.id) * (b.id - center.id) + (b.value - center.value) * (b.value - center.value);
	//cout << d1 << " " << d2 << endl;
	return d1 > d2;
}

TEMPLATE
bool GEOMETRY::PointCmp(POINT& center, POINT a, POINT b) {
	//GEOMETRY::POINT center;
	//center.value = center1.value;
	//center.id = center1.id;
	//GEOMETRY::center1.value = 0;
	//GEOMETRY::center1.id = 0;
	//cout << center.id<<" "<< center.value << endl;
	if (a.id >= 0 && b.id < 0)
		return true;
	if (a.id == 0 && b.id == 0)
		return a.value > b.value;
	/////////////////////////// vector product of OA and OB
	int det = (a.id - center.id) * (b.value - center.value) - (b.id - center.id) * (a.value - center.value);
	if (det < 0)
		return true;
	if (det > 0)
		return false;
	/////////////////////////////////// OA and OB in the same line, according to the distance to decide priority.
	int d1 = (a.id - center.id) * (a.id - center.id) + (a.value - center.value) * (a.value - center.value);
	int d2 = (b.id - center.id) * (b.id - center.value) + (b.value - center.value) * (b.value - center.value);
	return d1 > d2;
}

TEMPLATE
struct  GEOMETRY::CMP {
public:

	GEOMETRY::POINT center = NULL;
	GEOMETRY::CMP(GEOMETRY::POINT center_argu) : center(center_argu) {}

	bool operator() (const POINT& a, const POINT& b) {
		assert(center.id != NULL && center.value != NULL);
		return typename GEOMETRY::PointCmpCenter(center, a, b);
	}
};

TEMPLATE
bool GEOMETRY::Cmp(const POINT& a, const POINT& b) {
	//POINT center = center1;
	return PointCmp(a, b);
}

//************************************
// Method:Gravity
// Qualifier: Compute gravity center
// date:190307
// author:
//************************************
TEMPLATE
typename GEOMETRY::POINT GEOMETRY::Gravity(const vector<GEOMETRY::POINT>& const p) {
	assert(p.size() > 0);
	//ACCUMULATE_ID accumulate_id;
	//ACCUMULATE_VALUE accumulate_value;
	long double area = 0;
	int n = p.size();
	POINT center(0, 0);
	center.id = 0;
	center.value = 0;

	for (int i = 0; i < n - 1; i++) {
		area += (p[i].id * p[i + 1].value - p[i + 1].id * p[i].value) / 2;
		//cout << "area+: "<< area <<endl;
		center.id += (p[i].id * p[i + 1].value - p[i + 1].id * p[i].value) * (p[i].id + p[i + 1].id);
		center.value += (p[i].id * p[i + 1].value - p[i + 1].id * p[i].value) * (p[i].value + p[i + 1].value);
	}
	//cout << "area: " << area << endl;
	//cout << "add: "<< (p[n - 1].id*p[0].value - p[0].id*p[n - 1].value) / 2 <<endl;
	area += (p[n - 1].id * p[0].value - p[0].id * p[n - 1].value) / 2;

	if (area == 0 || area < (std::numeric_limits<double>::min)()) return getAverage(p, center); //190307 if area==0, points on a line.

	//cout << "area: "<< area <<endl;
	center.id += (p[n - 1].id * p[0].value - p[0].id * p[n - 1].value) * (p[n - 1].id + p[0].id);
	center.value += (p[n - 1].id * p[0].value - p[0].id * p[n - 1].value) * (p[n - 1].value + p[0].value);

	center.id /= 6.0 * area;
	center.value /= 6.0 * area;
	/*std::vector<POINT>::iterator it = p.begin();
	double sum_x = 0, sum_y = 0;
	for (auto&&i : p) {
		sum_x += i.x;
		sum_y += i.y;
	}*/
	//center.x = sum_x / p.size();
	//center.y = sum_y / p.size();
	//cout <<"************************" <<center.id << " " << center.value << endl;

	//auto[min, max] = minmax_element(p.begin(),p.end());

	return center;
}

//TEMPLATE
//void GEOMETRY::sort(vector <POINT>& vPoints) {
//
//}

TEMPLATE
void GEOMETRY::ClockwiseSortPoints(vector <GEOMETRY::POINT>& vPoints) {
	assert(vPoints.size() > 0);
	//GEOMETRY::center1 = Gravity(vPoints);
	//getPointsMaxDist(vPoints, Gravity(vPoints));//190311
	CMP Cmp(Gravity(vPoints));
	//Cmp.center = center1;

	/*for (auto&&i : vPoints) {
		cout << i.id << " " << i.value << "; ";
	}*/
	//cout << "center: [" << Cmp.center.id << ", " << Cmp.center.value << "]" << endl;
	sort(vPoints.begin(), vPoints.end(), Cmp);
}

TEMPLATE
int GEOMETRY::get_quadrant(const POINT& p) {
	int result = 4; //origin

	if (p.id > 0 && p.value > 0)
		return 1;
	else if (p.id < 0 && p.value > 0)
		return 2;
	else if (p.id < 0 && p.value < 0)
		return 3;
	//else 4th quadrant
	return result;
}

TEMPLATE
double GEOMETRY::get_clockwise_angle(const GEOMETRY::POINT& p) {
	double angle = 0.0;
	int quadrant = get_quadrant(p);

	/*making sure the quadrants are correct*/
	//cout << "Point: " << p << " is on the " << quadrant << " quadrant" << endl;

	/*calculate angle and return it*/
	angle = -atan2(p.id, -p.value);
	return angle;
}

TEMPLATE
bool GEOMETRY::compare_points(const POINT& a, const POINT& b) {
	/*if (a.x <= 0 && b.x > 0) return true;
	else if (a.x == 0 && b.x == 0) return a.y < b.y;
	else if (get_clockwise_angle(a) == get_clockwise_angle(b))
		return a.x*a.x + a.y*a.y > b.x*b.x + b.y*b.y;
	else*/
	//if (get_clockwise_angle(a) == get_clockwise_angle(b))
		//return a.x*a.x + a.y*a.y > b.x*b.x + b.y*b.y;
	return (get_clockwise_angle(a) < get_clockwise_angle(b));
}

//************************************
// Method:ComputePolygonArea
// Qualifier: Compute polygon area
// date:190225
// author:
//************************************
TEMPLATE
double GEOMETRY::ComputePolygonArea(const vector<POINT>& points) {//190225
	int point_num = points.size();
	if (point_num < 3) return 0.0;
	double s = points[0].value * (points[point_num - 1].id - points[1].id);
	for (int i = 1; i < point_num; ++i)
		s += points[i].value * (points[i - 1].id - points[(i + 1) % point_num].id);
	return fabs(s / 2.0);
}

TEMPLATE
struct GEOMETRY::COMPARE_POINTS {
	bool operator() (POINT a, POINT b) {
		return (GEOMETRY::get_clockwise_angle(a) < GEOMETRY::get_clockwise_angle(b));
	}
};

TEMPLATE
typename GEOMETRY::POINT GEOMETRY::nextToTop(stack<POINT>& S) {
	POINT p = S.top();
	S.pop();
	POINT res = S.top();
	S.push(p);
	return res;
}

TEMPLATE
void GEOMETRY::swap(POINT& p1, POINT& p2) {
	POINT temp = p1;
	p1 = p2;
	p2 = temp;
}

TEMPLATE
double GEOMETRY::distSq(POINT p1, POINT p2) {
	return (p1.id - p2.id) * (p1.id - p2.id) + (p1.value - p2.value) * (p1.value - p2.value);
}

TEMPLATE
int GEOMETRY::orientation(POINT p, POINT q, POINT r) {
	double val = (q.value - p.value) * (r.id - q.id) - (q.id - p.id) * (r.value - q.value);
	if (val == 0) return 0;  // colinear
	return (val > 0) ? 1 : 2; // clock:1 or counterclock:2 wise
}

TEMPLATE
int GEOMETRY::compare(const void* vp1, const void* vp2) {
	POINT* p1 = (POINT*)vp1;
	POINT* p2 = (POINT*)vp2;

	// Find orientation
	int o = orientation(p0, *p1, *p2);
	if (o == 0)
		return (distSq(p0, *p2) >= distSq(p0, *p1)) ? -1 : 1;

	return (o == 2) ? -1 : 1;
}

TEMPLATE
struct GEOMETRY::COMPARE_POLOR_POINTS {
public:

	GEOMETRY::POINT min_point;

	//COMPARE_POLOR_POINTS(double x, double y) : id(x), value(y) {}
	//COMPARE_POLOR_POINTS() : id(0), value(0) {}
	~COMPARE_POLOR_POINTS() {
		//id = NULL;
		//value = NULL;
	}

	bool operator() (POINT vp1, POINT vp2) {
		POINT* p1 = &vp1;
		POINT* p2 = &vp2;

		// Find orientation
		int o = GEOMETRY::orientation(min_point, *p1, *p2);
		if (o == 0)
			return (GEOMETRY::distSq(min_point, *p2) >= GEOMETRY::distSq(min_point, *p1)) ? true : false;

		return (o == 2) ? true : false;// 1 clockwise, 2 counterclockwise
	}

	/*int operator() (const void *vp1, const void *vp2) {
		POINT *p1 = (POINT *)vp1;
		POINT *p2 = (POINT *)vp2;

		 Find orientation
		int o = orientation(p0, *p1, *p2);
		if (o == 0)
			return (distSq(p0, *p2) >= distSq(p0, *p1)) ? -1 : 1;

		return (o == 2) ? -1 : 1;
	}
	int a;*/
};

//************************************
// Method:convexHull
// Qualifier: Get convex hull of polygon
// date:190301
// author:
//************************************
TEMPLATE
void GEOMETRY::convexHull(POINT points[], int n, std::vector <GEOMETRY::POINT>& points_convex) {
	// Find the bottommost point
	COMPARE_POLOR_POINTS compare_polor_points;
	double min_value = points[0].value;
	int min_id = 0;
	for (int i = 1; i < n; i++)
	{
		double value = points[i].value;

		// Pick the bottom-most or chose the left most point in case of tie
		if ((value < min_value) || (min_value == value && points[i].id < points[min_id].id))
			min_value = points[i].value, min_id = i;
	}

	// Place the bottom-most point at first position
	swap(points[0], points[min_id]);

	// Sort n-1 points with respect to the first point. A point p1 comes before p2 in sorted ouput if p2 has larger polar angle (in counterclockwise direction) than p1
	GEOMETRY::p0 = points[0];

	//vector<GEOMETRY::POINT> points_vector;

	//cout << "******  " << endl;
	//for (int i = 0; i < n; i++) {
	//	cout << points[i].id << " "<< points[i].value<<", ";
	//	points_vector.push_back(points[i]);
	//}
	//cout << endl;
	compare_polor_points.min_point = points[0];
	//qsort(&points[1], n - 1, sizeof(POINT), compare_polor_points);
	//sort(points_vector.begin()+1, points_vector.end(), compare_polor_points);

	sort(points + 1, points + n, compare_polor_points);

	/*for (auto&&i: points_vector) {
		cout << i.id << " " << i.value << ", ";
	}

	for (int i = 0; i < n; i++) {
		assert(points_vector.at(i).id == points[i].id && points_vector.at(i).value == points[i].value);
	}*/

	//cout << endl;
	// If two or more points make same angle with p0, Remove all but the one that is farthest from p0.
	// Remember that, in above sorting, our criteria was to keep the farthest point at the end when more than one points have same angle.
	int m = 1; // Initialize size of modified array
	for (int i = 1; i < n; i++) {
		// Keep removing i while angle of i and i+1 is same with respect to p0
		while (i < n - 1 && orientation(p0, points[i], points[i + 1]) == 0) {
			//cout << "coline" << endl;
			i++;
		}

		points[m] = points[i];
		m++;  // Update size of modified array
	}

	// If modified array of points has less than 3 points,
	// convex hull is not possible
	if (m < 3) return;

	// Create an empty stack and push first three points
	// to it.
	stack<GEOMETRY::POINT> S;
	S.push(points[0]);
	S.push(points[1]);
	S.push(points[2]);

	// Process remaining n-3 points
	for (int i = 3; i < m; i++)
	{
		// Keep removing top while the angle formed by
		// points next-to-top, top, and points[i] makes
		// a non-left turn
		while (orientation(nextToTop(S), S.top(), points[i]) != 2)
			S.pop();
		S.push(points[i]);
	}

	// Now stack has the output points, print contents of stack

	assert(!S.empty());

	//cout << "points_convex: ";
	while (!S.empty()) {
		GEOMETRY::POINT p = S.top();
		points_convex.push_back(S.top());
		//cout << "(" << p.id << ", " << p.value << "), ";
		S.pop();
	}

	//cout << endl;
}

//************************************
// Method:isHollow
// Qualifier: Judge hollow
// date:190306
// author:
//************************************
TEMPLATE
bool GEOMETRY::isHollow(std::vector<GEOMETRY::POINT>& curveloopPoints) {
	auto num = curveloopPoints.size();
	double angleSum = 0.0;
	for (int i = 0; i < num; i++) {
		double e1;
		if (i == 0) {
			e1 = curveloopPoints[num - 1].value - curveloopPoints[i].value;
		}
		else {
			e1 = curveloopPoints[i - 1].value - curveloopPoints[i].value;
		}
		double e2;
		if (i == num - 1) {
			e2 = curveloopPoints[0].value - curveloopPoints[i].value;
		}
		else {
			e2 = curveloopPoints[i + 1].value - curveloopPoints[i].value;
		}
		//标准化并计算点乘
		//e1.normalize(); e2.normalize();
		//double mdot = e1 % e2;
		double mdot = std::fmod(e1, e2);
		//计算较小值
		double theta = acos(mdot);
		//加和
		angleSum += theta;
	}
	//计算内角和
	double convexAngleSum = double((num - 2)) * 3.1415926;
	//判断凹凸性
	if (angleSum < (convexAngleSum - (double(num) * 0.00001)))
	{/*
		if (HollowPoints.size()>0)
		{*/
		//m_IsHollow = true;
	//}
		return true;//是凹
	}
	return false;//otherwise is convex
}

//************************************
// Method:isConvex
// Qualifier: Judge convex
// date:190306
// author:
//************************************
//TEMPLATE
//int GEOMETRY::isConvex(POINT *p, int n){//n为点数
//
//	int i, j, k;
//	int flag = 0;
//	double z;
//	if (n < 3)
//		return 0;
//
//	for (i = 0; i < n; i++) {
//		j = (i + 1) % n;
//		k = (i + 2) % n;
//		z = (p[j].id - p[i].id) * (p[k].value - p[j].value);
//		z -= (p[j].value - p[i].value) * (p[k].id - p[j].id);
//
//		if (z < 0)
//			flag |= 1;
//		else if (z > 0)
//			flag |= 2;
//		if (flag == 3)
//			return -1;//凸集
//	}
//	if (flag != 0) //凹集
//		return 1;
//	else
//		return 0; //曲线不合法
//}

//************************************
// Method:isConvex
// Qualifier: Judge convex
// date:190306
// author:
//************************************
TEMPLATE
int GEOMETRY::isConvex(std::vector<GEOMETRY::POINT>& p) {
	int i, j, k;
	int n = p.size();
	int flag = 0;
	double z;
	if (n < 3)
		return 0;

	for (i = 0; i < n; i++) {
		j = (i + 1) % n;
		k = (i + 2) % n;
		z = (p[j].id - p[i].id) * (p[k].value - p[j].value);
		z -= (p[j].value - p[i].value) * (p[k].id - p[j].id);

		if (z < 0)
			flag |= 1;
		else if (z > 0)
			flag |= 2;
		if (flag == 3)
			return -1;//凸集
	}
	if (flag != 0) //凹集
		return 1;
	else
		return 0; //曲线不合法
}

TEMPLATE
template<typename T>
inline T GEOMETRY::area_of_a_circle(T r) {
	//using boost::math::constants::pi;
	//return pi<T>() * r * r;
}

//************************************
// Method:get2PointsDistSquare
// Qualifier:Get the square distance between 2 points.
// date:190308
// author:
//************************************
TEMPLATE
double GEOMETRY::get2PointsDist(const GEOMETRY::POINT& const point_1, const GEOMETRY::POINT& const point_2) {//190308
	//3 options for distance function

	typedef boost::geometry::model::d2::point_xy<double> point_type;
	point_type p1(point_1.id, point_1.value);
	point_type p2(point_2.id, point_2.value);

	auto dist1 = sqrt((point_1.id - point_2.id) * (point_1.id - point_2.id) + (point_1.value - point_2.value) * (point_1.value - point_2.value));//distance1

	auto dist2 = std::hypot(point_1.id - point_2.id, point_1.value - point_2.value);//distance2

	auto dist3 = boost::geometry::distance(p1, p2);//distance3

	assert(dist1 == dist2 && dist2 == dist3 && dist1 == dist3);

	return dist3;
}

//************************************
// Method:getPointsMaxDist
// Qualifier:Get the maximum distance between center poit and other points.
// date:190308
// author:
//************************************
TEMPLATE
double GEOMETRY::getPointsMaxDist(const std::vector<GEOMETRY::POINT>& const points, const  GEOMETRY::POINT& const center_point) {//190308
	double max_dist = DBL_MIN;

	for (auto&& i : points) {
		assert(get2PointsDist(center_point, i) == std::hypot(i.id - center_point.id, i.value - center_point.value));
		max_dist = max_dist < get2PointsDist(center_point, i) ? get2PointsDist(center_point, i) : max_dist;
	}

	return max_dist;
}

//************************************
// Method:sgn(double x)
// Qualifier:
// date:190311
// author:
//************************************
TEMPLATE
int GEOMETRY::sgn(double x) {
	if (fabs(x) < EPS)
		return 0;
	return x < 0 ? -1 : 1;
}

//************************************
// Method:get_distance
// Qualifier:
// date:190311
// author:
//************************************
TEMPLATE
double GEOMETRY::get_distance(const POINT a, const POINT b)//两点之间的距离
{
	return sqrt((a.id - b.id) * (a.id - b.id) + (a.value - b.value) * (a.value - b.value));
}

//************************************
// Method:get_circle_center
// Qualifier:
// date:190311
// author:
//************************************
TEMPLATE
typename GEOMETRY::POINT GEOMETRY::get_circle_center(const GEOMETRY::POINT a, const GEOMETRY::POINT b, const GEOMETRY::POINT c) {//得到triangle外接圆的圆心
	GEOMETRY::POINT center;
	double a1 = b.id - a.id;
	double b1 = b.value - a.value;
	double c1 = (a1 * a1 + b1 * b1) / 2.0;
	double a2 = c.id - a.id;
	double b2 = c.value - a.value;
	double c2 = (a2 * a2 + b2 * b2) / 2.0;
	double d = a1 * b2 - a2 * b1;
	center.id = a.id + (c1 * b2 - c2 * b1) / d;
	center.value = a.value + (a1 * c2 - a2 * c1) / d;
	return center;
}

//************************************
// Method:get_circle_center
// Qualifier:
// date:190311
// author:
//************************************
TEMPLATE
//p表示定点, n表示顶点的个数, c代表最小覆盖圆圆心, r是半径
void GEOMETRY::min_cover_circle(const std::vector<GEOMETRY::POINT>& const points, GEOMETRY::POINT& const center_point, double& const circle_radius)//找最小覆盖圆(这里没有用全局变量p[], 因为是为了封装一个函数便于调用)
{
	/*typedef boost::geometry::model::d2::point_xy<double> point_type;
	point_type point_boost;
	point_type center_boost;*/
	//static std::default_random_engine randGen((std::random_device())());
	//std::shuffle(points.begin(), points.end(), randGen);
	//shuffle(points.begin(), points.end(),seed);//随机函数,使用了之后使程序更快点,也可以不用

	center_point = points[0];
	circle_radius = 0;
	for (int i = 1; i < points.size(); i++) {
		/*boost::geometry::assign_values(point_boost, points[i].id, points[i].value);
		boost::geometry::assign_values(center_boost, center_point.id, center_point.value);
		boost::geometry::distance(point_boost, center_boost);*/
		//boost::geometry::distance(boost::geometry::assign_values(point_boost, points[i]), boost::geometry::assign_values(center_boost, center_point));

		if (sgn(get2PointsDist(points[i], center_point) - circle_radius) > 0)//如果p[i]在当前圆的外面, 那么以当前点为圆心开始找
		{
			center_point = points[i];//圆心为当前点
			circle_radius = 0;//这时候这个圆只包括他自己.所以半径为0
			for (int j = 0; j < i; j++)//找它之前的所有点
			{
				if (sgn(get2PointsDist(points[j], center_point) - circle_radius) > 0)//如果之前的点有不满足的, 那么就是以这两点为直径的圆
				{
					center_point.id = (points[i].id + points[j].id) / 2.0;
					center_point.value = (points[i].value + points[j].value) / 2.0;
					circle_radius = get2PointsDist(points[j], center_point);
					for (int k = 0; k < j; k++)
					{
						if (sgn(get2PointsDist(points[k], center_point) - circle_radius) > 0)//找新作出来的圆之前的点是否还有不满足的, 如果不满足一定就是三个点都在圆上了
						{
							vector<GEOMETRY::POINT> points_on_circle = { points[i], points[j], points[k] };

							/*for (auto&&i : points_on_circle) {
								cout << i.id << " " << i.value << "; ";
							}
							cout << endl;*/

							//auto temp_center = Gravity(points_on_circle);
							center_point = get_circle_center(points[i], points[j], points[k]);

							//cout << "Gravity: " << temp_center.id << " " << temp_center.value << ";    center: " << center_point.id << " " << center_point.value <<endl;

							//auto temp_radius = get2PointsDist(points[i], temp_center);
							circle_radius = get2PointsDist(points[i], center_point);
							//cout << "Gravity radius: "<< temp_radius << "; center radius:"<<  circle_radius <<endl;
							//assert(center_point.id == temp_center.id && center_point.value == temp_center.value);
						}
					}
				}
			}
		}
	}
}



#endif

