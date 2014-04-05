#ifndef _GCM_AREA_H
#define _GCM_AREA_H 1

#include "../../datatypes/ElasticNode.h"
#include "../GCMException.h"

class Area
{
public:
	virtual bool isInArea( ElasticNode* cur_node ) = 0;
};

#endif
