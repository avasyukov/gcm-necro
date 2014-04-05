#ifndef _GCM_VOLUME_CALCULATOR_H
#define _GCM_VOLUME_CALCULATOR_H  1

#include <vector>
using std::vector;

#include "../datatypes/ElasticNode.h"
#include "../datatypes/ElasticMatrix3D.h"

class VolumeCalculator
{
public:
	VolumeCalculator();
	~VolumeCalculator();
	virtual void do_calc(ElasticNode* new_node, ElasticMatrix3D* matrix, float* values[]) = 0;

protected:

};

#endif
