#include "SPHConnector.h"

SPHConnector::SPHConnector() { };

SPHConnector::~SPHConnector() { };

// IN - coords - x, y, z - absolute values
// OUT - values - vx, vy, vz, sxx, sxy, sxz, syy, syz, szz (we do not need lambda, mu, rho, yield_limit for these 'virtual sph nodes')
void SPHConnector::interpolate( float* coords, float* values )
{
	values[2] = 20;
};
