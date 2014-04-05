#include "ContactCondition.h"

void ContactCondition::do_calc(float time, ElasticNode* cur_node, ElasticNode* virt_node, ElasticNode* new_node, ElasticMatrix3D* matrix, float* values[], bool inner[], ElasticMatrix3D* virt_matrix, float* virt_values[], bool virt_inner[], float outer_normal[])
{
	calc->do_calc(cur_node, virt_node, new_node, matrix, values, inner, virt_matrix, virt_values, virt_inner, outer_normal, form->calcMagnitudeNorm(time, cur_node->coords, area) );
};
