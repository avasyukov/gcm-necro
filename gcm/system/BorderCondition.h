#ifndef _GCM_BORDER_CONDITION_H
#define _GCM_BORDER_CONDITION_H 1

#include "./areas/Area.h"
#include "./forms/PulseForm.h"
#include "../methods/BorderCalculator.h"

class BorderCondition
{
public:
	Area* area;
	PulseForm* form;
	BorderCalculator* calc;
	void do_calc(float time, float* cur_coords, ElasticNode* new_node, ElasticMatrix3D* matrix, float* values[], bool inner[], float outer_normal[]);
};

#endif
