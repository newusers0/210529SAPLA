#pragma once

#ifndef PCH_H
#define PCH_H

//#include "targetver.h"

#include "stdafx.h"
#include <stdio.h>
#include <tchar.h>
//#include <vld.h>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <stack>
#include <string>
#include <cassert>
#include <iomanip>
#include <functional>
#include <numeric>
#include <time.h>
#include <chrono>
#include <ctime>
#include <Windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <cmath>
#include <math.h>
#include <complex>
#include <random>
#include <type_traits>
#include <iosfwd>
#include <limits>

/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/polygon/polygon.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>//200804 print system time
namespace bg = boost::geometry;
namespace bp = boost::polygon;
using namespace boost::polygon::operators;
namespace gtl = boost::polygon;
typedef boost::geometry::model::d2::point_xy<double> boost_point_type;
/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
typedef double DataType;
typedef double ElementType;
/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
#define TEMPLATE template<typename DataType>
#define TOOL SHARE_TOOL<DataType>
#define GEOMETRY GEOMETRY_TOOL<DataType>
#define APCA_QUAL CAPCA<DataType>
#define APCA_KNN_QUAL CAPCA_KNN<DataType>
#define PLA_QUAL  CPLA<DataType>
#define CHEBYSHEV_QUAL  CCHEBYSHEV<DataType>
#define DFT CDFT<DataType>
#define APLA CAPLA<DataType>//181112
#define MULTI MULTI_DIMENSION<DataType>
#define Evaluation EVALUATION<DataType>

#define INF std::numeric_limits<double>::infinity()
#define EPS 1e-8//190311 for min circle cover
#define MIN_D 1e-4
/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

#include "./lib/RTree.h"
typedef RTree<DataType, ElementType> RTREE;

static int PAA_or_APCA = DBL_MAX;


// TODO: add headers that you want to pre-compile here

#endif //PCH_H