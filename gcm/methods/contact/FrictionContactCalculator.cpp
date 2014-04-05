#include "FrictionContactCalculator.h"
#include "AdhesionContactCalculator.h"
#include "../border/ExternalForceCalculator.h"
#include "../border/FreeBorderCalculator.h"
#include "../../datatypes/ElasticNode.h"

FrictionContactCalculator::FrictionContactCalculator()
{
	U_gsl = gsl_matrix_alloc (18, 18);
	om_gsl = gsl_vector_alloc (18);
	x_gsl = gsl_vector_alloc (18);
	p_gsl = gsl_permutation_alloc (18);
};

FrictionContactCalculator::~FrictionContactCalculator()
{
	gsl_matrix_free(U_gsl);
	gsl_vector_free(om_gsl);
	gsl_vector_free(x_gsl);
	gsl_permutation_free(p_gsl);
};

void FrictionContactCalculator::do_calc(ElasticNode* cur_node,ElasticNode* virt_node, ElasticNode* new_node, ElasticMatrix3D* matrix, float* values[], bool inner[], ElasticMatrix3D* virt_matrix, float* virt_values[], bool virt_inner[], float outer_normal[], float scale)
{
	float local_n[3][3];
	local_n[0][0] = outer_normal[0];
	local_n[0][1] = outer_normal[1];
	local_n[0][2] = outer_normal[2];

	qm.create_local_basis(local_n[0], local_n[1], local_n[2]);

	//---------------------------------------Check if nodes fall apart
	float vel_rel[3] = {cur_node->values[0]-virt_node->values[0],
		   	    cur_node->values[1]-virt_node->values[1],
			    cur_node->values[2]-virt_node->values[2]};
	float force_cur[3] = {cur_node->values[3]*outer_normal[0]+cur_node->values[4]*outer_normal[1]+cur_node->values[5]*outer_normal[2],
			  cur_node->values[4]*outer_normal[0]+cur_node->values[6]*outer_normal[1]+cur_node->values[7]*outer_normal[2],
			  cur_node->values[5]*outer_normal[0]+cur_node->values[7]*outer_normal[1]+cur_node->values[8]*outer_normal[2]};
	float   dsxx =  -virt_node->values[3] + cur_node->values[3],
		dsxy =  -virt_node->values[4] + cur_node->values[4],
		dsxz =  -virt_node->values[5] + cur_node->values[5],
		dsyy =  -virt_node->values[6] + cur_node->values[6],
		dsyz =  -virt_node->values[7] + cur_node->values[7],
		dszz =  -virt_node->values[8] + cur_node->values[8];
	float rho = cur_node->rho,
		c1 = sqrt( (cur_node->la + 2 * cur_node->mu) / rho ),
		c2 = sqrt( cur_node->mu / rho );
	float vel_norm_abs =  qm.scalar_product(vel_rel[0],vel_rel[1],vel_rel[2],outer_normal[0],outer_normal[1],outer_normal[2]);
	float vel_norm_sign = 1.0;
	if (vel_norm_sign < 0.0) vel_norm_sign = -1.0;
	float vel_tang[3] = {vel_rel[0] - vel_norm_abs*outer_normal[0],
				vel_rel[1] - vel_norm_abs*outer_normal[1],
				vel_rel[2] - vel_norm_abs*outer_normal[2]};
	float vel_tang_abs = sqrt(qm.scalar_product(vel_tang[0],vel_tang[1],vel_tang[2], vel_tang[0],vel_tang[1],vel_tang[2]));
	float force_cur_abs = qm.scalar_product(force_cur[0],force_cur[1],force_cur[2],outer_normal[0],outer_normal[1],outer_normal[2]);
	float k = 0.1;
	
	float eps = 0.000000005;
	bool free_border = false, frict = false, adh_frict = false, sli_frict = false;

	if (vel_norm_abs > eps) 			//first check relative speed
		free_border = true;
	else if (vel_norm_abs > - eps)			//if relative speed is small, we check force
	{
		if (force_cur_abs < eps)
			free_border = true;
		else frict = true;
	}
	else frict = true;

	//if (adh_frict)
        //{  
//                AdhesionContactCalculator* acc = new AdhesionContactCalculator();
  //              acc->do_calc(cur_node,new_node,virt_node,matrix,values,inner,virt_matrix,virt_values,virt_inner,outer_normal,scale);
    //            delete acc;
      //          return;
        //}

	if (false)//free_border) 		
	{
		FreeBorderCalculator *fbc = new FreeBorderCalculator();
		fbc->do_calc(new_node, matrix, values, inner, outer_normal, scale);
		delete fbc;
//		LOG_INFO("Friction uses free border");
		return;
	}
	
      
		ElasticNode* nd = new ElasticNode();
                AdhesionContactCalculator* acc = new AdhesionContactCalculator();
                acc->do_calc(cur_node,virt_node,nd,matrix,values,inner,virt_matrix,virt_values,virt_inner,outer_normal,scale);
                delete acc;
	        float frc_nd[3] = 
			{	nd->values[3]*outer_normal[0]+nd->values[4]*outer_normal[1]+nd->values[5]*outer_normal[2],
                          	nd->values[4]*outer_normal[0]+nd->values[6]*outer_normal[1]+nd->values[7]*outer_normal[2],
                          	nd->values[5]*outer_normal[0]+nd->values[7]*outer_normal[1]+nd->values[8]*outer_normal[2]};
		float frc_tan = qm.scalar_product(vel_tang[0],vel_tang[1],vel_tang[2],frc_nd[0],frc_nd[1],frc_nd[2])/vel_tang_abs;
		float frc_nrm = qm.scalar_product(outer_normal[0],outer_normal[1],outer_normal[2],frc_nd[0],frc_nd[1],frc_nd[2]);
		if (k*frc_nrm > frc_tan)
			sli_frict = true;
		else
			sli_frict = false;
    

	if (sli_frict)
        {
                ExternalForceCalculator *efc = new ExternalForceCalculator();
		float tp[3] = {dsxx*outer_normal[0]+dsxy*outer_normal[1]+dsxz*outer_normal[2],
				dsxy*outer_normal[0]+dsyy*outer_normal[1]+dsyz*outer_normal[2],
				dsxz*outer_normal[0]+dsyz*outer_normal[1]+dszz*outer_normal[2]};
		float sn =  (- 0.5*rho*c1*vel_norm_abs - 0.5*qm.scalar_product(tp[0],tp[1],tp[2],outer_normal[0],outer_normal[1],outer_normal[2]));
		//	0.5*(rho*c1*vel_norm_abs);
		if (sn > 0.0) sn = -sn;
		efc->set_parameters(sn,k*sn,vel_tang[0]/vel_tang_abs,vel_tang[1]/vel_tang_abs,vel_tang[2]/vel_tang_abs);
                efc->do_calc(new_node, matrix, values, inner, outer_normal, scale);
                delete efc;
                return;
        }

        else //if (adh_frict)
        {
                AdhesionContactCalculator* acc = new AdhesionContactCalculator();
                acc->do_calc(cur_node,new_node,virt_node,matrix,values,inner,virt_matrix,virt_values,virt_inner,outer_normal,scale);
                delete acc;
                return;
        }

	//--------------------------------------------------------------------



};
