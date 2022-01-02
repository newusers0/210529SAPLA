#pragma once

#pragma once
#ifndef CAPLA_PUB_H
#define CAPLA_PUB_H

#include "pch.h"
#include "SHARE_TOOL.h"
#include "GEOMETRY_TOOL.h"
#include "lib/doublyLinkedList.h"
#include "lib/SmallestEnclosingCircle.hpp"
#include "CPLA.h"
//#include "./lib/RTree.h"//210618

class CAPLA_PUB : public RTREE, virtual public GEOMETRY, virtual public TOOL, virtual public PLA_QUAL{


	//210509 get triangle area by line A1A2 and B1B2, they have same length
	template<typename T>
	inline long double get_triangle_area_pub(const T& const segment_1, const T& const segment_2);

	//***************************************************************
	// Method:get_distance_SAPLA_pub
	// Qualifier:  Lower bound distance between SAPLA (Approximation)
	// Input:
	// Output:
	// notice:
	// date:210506
	// author:
	//***************************************************************
	template<typename T, typename Y>
	long double get_distance_SAPLA_pub(const vector<T>& const original_time_series_vector_1, const vector<T>& const original_time_series_vector_2, const DoublyLinkedList<Y>& const doubly_linked_list_1, const DoublyLinkedList<Y>& const doubly_linked_list_2);

};

#endif

