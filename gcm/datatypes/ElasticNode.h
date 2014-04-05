#ifndef _GCM_ELASTIC_NODE_H
#define _GCM_ELASTIC_NODE_H  1

#include <vector>
using std::vector;

#include "Node.h"

class TetrMesh_1stOrder;
class VolumeCalculator;
class BorderCondition;
class ContactCondition;

class ElasticNode : public Node
{
/* Inherited from Node 
*	int zone_num;
*	int local_num;
*	int remote_num;
*	int absolute_num;
*	int placement_type;
*	float coords[3];
*	float fixed_coords[3];
*/
public:
	union
	{
		float values[13];
		struct
		{
			float vx;
			float vy;
			float vz;
			float sxx;
			float sxy;
			float sxz;
			float syy;
			float syz;
			float szz;
			float la;	// TODO If la and mu should be replaced with some ID of rheology class???
			float mu;	//
			float rho;
			float yield_limit;
		};
	};

	union
	{
		float destruction_criterias[8];
		struct
		{
			float max_compression;
			float max_tension;
			float max_shear;
			float max_deviator;
			float max_compression_history;
			float max_tension_history;
			float max_shear_history;
			float max_deviator_history;
		};
	};

	vector<int>* elements;
	vector<int>* border_elements;
	TetrMesh_1stOrder* mesh;

	VolumeCalculator* volume_calculator;
	BorderCondition* border_condition;
	ContactCondition* contact_condition;
protected:
	// TODO should we switch from vector to memmory block with offsets in it?
	// int elems_offset;
	// int elems_size;
	// TODO How to deal with temporary variables???
	// int element_for_interpolation;
	// float random_axis[3];
	// float maxL;
	// .. and so on
};

#endif
