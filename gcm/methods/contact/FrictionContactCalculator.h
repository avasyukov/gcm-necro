#ifndef _GCM_CONTACT_FRICTION_CALCULATOR_H
#define _GCM_CONTACT_FRICTION_CALCULATOR_H  1

#include "../ContactCalculator.h"
#include "../../system/quick_math.h"
#include <gsl/gsl_linalg.h>

class FrictionContactCalculator : public ContactCalculator
{
public:
	FrictionContactCalculator();
	~FrictionContactCalculator();
	void do_calc(ElasticNode* cur_node,ElasticNode* vitr_node, ElasticNode* new_node, ElasticMatrix3D* matrix, float* values[], bool inner[], ElasticMatrix3D* virt_matrix, float* virt_values[], bool virt_inner[], float outer_normal[], float scale);

protected:

private:
	quick_math qm;
	// Used for border calculation
	gsl_matrix *U_gsl;
	gsl_vector *om_gsl;
	gsl_vector *x_gsl;
	gsl_permutation *p_gsl;
};

#endif
