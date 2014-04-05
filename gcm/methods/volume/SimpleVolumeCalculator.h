#ifndef _GCM_SIMPLE_VOLUME_CALCULATOR_H
#define _GCM_SIMPLE_VOLUME_CALCULATOR_H  1

#include "../VolumeCalculator.h"

class SimpleVolumeCalculator : public VolumeCalculator
{
public:
	SimpleVolumeCalculator();
	~SimpleVolumeCalculator();
	void do_calc(ElasticNode* new_node, ElasticMatrix3D* matrix, float* values[]);

protected:

};

#endif
