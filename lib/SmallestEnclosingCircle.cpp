

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include "SmallestEnclosingCircle.hpp"

using std::size_t;
using std::vector;
using std::max;
using std::min;


/*---- Members of struct Point ----*/

Point_Enclose Point_Enclose::subtract(const Point_Enclose &p) const {
	return Point_Enclose{x - p.x, y - p.y};
}


double Point_Enclose::distance(const Point_Enclose &p) const {
	return std::hypot(x - p.x, y - p.y);
}


double Point_Enclose::cross(const Point_Enclose &p) const {
	return x * p.y - y * p.x;
}


/*---- Members of struct Circle ----*/

const Circle Circle::INVALID{Point_Enclose{0, 0}, -1};

const double Circle::MULTIPLICATIVE_EPSILON = 1 + 1e-14;


bool Circle::contains(const Point_Enclose &p) const {
	return c.distance(p) <= r * MULTIPLICATIVE_EPSILON;
}


bool Circle::contains(const vector<Point_Enclose> &ps) const {
	for (const Point_Enclose &p : ps) {
		if (!contains(p))
			return false;
	}
	return true;
}


/*---- Smallest enclosing circle algorithm ----*/

static Circle makeSmallestEnclosingCircleOnePoint (const vector<Point_Enclose> &points, size_t end, const Point_Enclose &p);
static Circle makeSmallestEnclosingCircleTwoPoints(const vector<Point_Enclose> &points, size_t end, const Point_Enclose &p, const Point_Enclose &q);

static std::default_random_engine randGen((std::random_device())());


// Initially: No boundary points known
Circle makeSmallestEnclosingCircle(const vector<Point_Enclose> &points) {
	// Clone list to preserve the caller's data, randomize order
	vector<Point_Enclose> shuffled = points;
	std::shuffle(shuffled.begin(), shuffled.end(), randGen);
	
	// Progressively add points to circle or recompute circle
	Circle c = Circle::INVALID;
	for (size_t i = 0; i < shuffled.size(); i++) {
		const Point_Enclose &p = shuffled.at(i);
		if (c.r < 0 || !c.contains(p))
			c = makeSmallestEnclosingCircleOnePoint(shuffled, i + 1, p);
	}
	return c;
}


// One boundary point known
static Circle makeSmallestEnclosingCircleOnePoint(const vector<Point_Enclose> &points, size_t end, const Point_Enclose &p) {
	Circle c{p, 0};
	for (size_t i = 0; i < end; i++) {
		const Point_Enclose &q = points.at(i);
		if (!c.contains(q)) {
			if (c.r == 0)
				c = makeDiameter(p, q);
			else
				c = makeSmallestEnclosingCircleTwoPoints(points, i + 1, p, q);
		}
	}
	return c;
}


// Two boundary points known
static Circle makeSmallestEnclosingCircleTwoPoints(const vector<Point_Enclose> &points, size_t end, const Point_Enclose &p, const Point_Enclose &q) {
	Circle circ = makeDiameter(p, q);
	Circle left  = Circle::INVALID;
	Circle right = Circle::INVALID;
	
	// For each point not in the two-point circle
	Point_Enclose pq = q.subtract(p);
	for (size_t i = 0; i < end; i++) {
		const Point_Enclose &r = points.at(i);
		if (circ.contains(r))
			continue;
		
		// Form a circumcircle and classify it on left or right side
		double cross = pq.cross(r.subtract(p));
		Circle c = makeCircumcircle(p, q, r);
		if (c.r < 0)
			continue;
		else if (cross > 0 && (left.r < 0 || pq.cross(c.c.subtract(p)) > pq.cross(left.c.subtract(p))))
			left = c;
		else if (cross < 0 && (right.r < 0 || pq.cross(c.c.subtract(p)) < pq.cross(right.c.subtract(p))))
			right = c;
	}
	
	// Select which circle to return
	if (left.r < 0 && right.r < 0)
		return circ;
	else if (left.r < 0)
		return right;
	else if (right.r < 0)
		return left;
	else
		return left.r <= right.r ? left : right;
}


Circle makeDiameter(const Point_Enclose &a, const Point_Enclose &b) {
	Point_Enclose c{(a.x + b.x) / 2, (a.y + b.y) / 2};
	return Circle{c, max(c.distance(a), c.distance(b))};
}


Circle makeCircumcircle(const Point_Enclose &a, const Point_Enclose &b, const Point_Enclose &c) {
	// Mathematical algorithm from Wikipedia: Circumscribed circle
	double ox = (min(min(a.x, b.x), c.x) + max(min(a.x, b.x), c.x)) / 2;
	double oy = (min(min(a.y, b.y), c.y) + max(min(a.y, b.y), c.y)) / 2;
	double ax = a.x - ox,  ay = a.y - oy;
	double bx = b.x - ox,  by = b.y - oy;
	double cx = c.x - ox,  cy = c.y - oy;
	double d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2;
	if (d == 0)
		return Circle::INVALID;
	double x = ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + (cx*cx + cy*cy) * (ay - by)) / d;
	double y = ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + (cx*cx + cy*cy) * (bx - ax)) / d;
	Point_Enclose p{ox + x, oy + y};
	double r = max(max(p.distance(a), p.distance(b)), p.distance(c));
	return Circle{p, r};
}
