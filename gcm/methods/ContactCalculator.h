#ifndef _GCM_CONTACT_CALCULATOR_H
#define _GCM_CONTACT_CALCULATOR_H  1

#include <vector>
using std::vector;

#include "../datatypes/ElasticNode.h"
#include "../datatypes/ElasticMatrix3D.h"

class ContactCalculator
{
public:
	ContactCalculator();
	~ContactCalculator();
	virtual void do_calc(ElasticNode* cur_node,ElasticNode* vitr_node, ElasticNode* new_node, ElasticMatrix3D* matrix, float* values[], bool inner[], ElasticMatrix3D* virt_matrix, float* virt_values[], bool virt_inner[], float outer_normal[], float scale) = 0;

protected:

};

#endif
