#pragma once

#ifndef EVALUATION_H
#define EVALUATION_H

#include "pch.h"
#include "SHARE_TOOL.h"
#include "APLA_ICDE07.h"
#include "CAPLA.h"
#include "MULTI_DIMENSION.h"

TEMPLATE
class EVALUATION : virtual public TOOL, virtual public APLA
{
public:
	void evaluate_multi_KNN_speed();
};

#include "EVALUATION.cpp"

#endif

