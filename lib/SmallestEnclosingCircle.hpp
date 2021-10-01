
#pragma once

#include <vector>


struct Point_Enclose final {
	
	public: double x;
	public: double y;
	
	
	public: Point_Enclose subtract(const Point_Enclose &p) const;
	
	public: double distance(const Point_Enclose &p) const;
	
	// Signed area / determinant thing
	public: double cross(const Point_Enclose &p) const;
	
};


struct Circle final {
	
	public: static const Circle INVALID;
	
	private: static const double MULTIPLICATIVE_EPSILON;
	
	
	public: Point_Enclose c;   // Center
	public: double r;  // Radius
	
	
	public: bool contains(const Point_Enclose &p) const;
	
	public: bool contains(const std::vector<Point_Enclose> &ps) const;
	
};



Circle makeSmallestEnclosingCircle(const std::vector<Point_Enclose> &points);


Circle makeDiameter(const Point_Enclose &a, const Point_Enclose &b);
Circle makeCircumcircle(const Point_Enclose &a, const Point_Enclose &b, const Point_Enclose &c);
