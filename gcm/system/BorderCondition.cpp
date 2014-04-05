#include "BorderCondition.h"

void BorderCondition::do_calc(float time, float* cur_coords, ElasticNode* new_node, ElasticMatrix3D* matrix, float* values[], bool inner[], float outer_normal[])
{
	calc->do_calc(new_node, matrix, values, inner, outer_normal, form->calcMagnitudeNorm(time, cur_coords, area) );
};
