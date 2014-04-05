#ifndef _GCM_BORDER_CALCULATOR_H
#define _GCM_BORDER_CALCULATOR_H  1

#include <vector>
using std::vector;

#include "../datatypes/ElasticNode.h"
#include "../datatypes/ElasticMatrix3D.h"

class BorderCalculator
{
public:
	BorderCalculator();
	~BorderCalculator();
	virtual void do_calc(ElasticNode* new_node, ElasticMatrix3D* matrix, float* values[], bool inner[], float outer_normal[], float scale) = 0;

protected:

};

#endif
