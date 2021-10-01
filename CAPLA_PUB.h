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


};

#endif

