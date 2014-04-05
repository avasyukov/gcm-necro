#ifndef _GCM_SPH_H
#define _GCM_SPH_H  1

#include "../system/LoggerUser.h"

class SPHConnector : protected LoggerUser
{
public:
	SPHConnector();
	~SPHConnector();
	void interpolate( float* coords, float* values );
};

#endif

