#ifndef _GCM_CONTACT_CONDITION_H
#define _GCM_CONTACT_CONDITION_H 1

#include "./areas/Area.h"
#include "./forms/PulseForm.h"
#include "../methods/ContactCalculator.h"

class ContactCondition
{
public:
	Area* area;
	PulseForm* form;
	ContactCalculator* calc;
	void do_calc(float time, ElasticNode* cur_node, ElasticNode* virt_node, ElasticNode* new_node, ElasticMatrix3D* matrix, float* values[], bool inner[], ElasticMatrix3D* virt_matrix, float* virt_values[], bool virt_inner[], float outer_normal[]);
};

#endif
