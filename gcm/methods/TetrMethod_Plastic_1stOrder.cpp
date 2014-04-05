#include "TetrMethod_Plastic_1stOrder.h"

GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis()
{
	num_method_type.assign("1st order interpolation on tetr mesh (with plasticity) - axis rotation");
	// TODO - do we really need 3 matrixes? May be we can use just one to store 'current stage'?
	elastic_matrix3d[0] = new ElasticMatrix3D();
	elastic_matrix3d[1] = new ElasticMatrix3D();
	elastic_matrix3d[2] = new ElasticMatrix3D();
	virt_elastic_matrix3d[0] = new ElasticMatrix3D();
	virt_elastic_matrix3d[1] = new ElasticMatrix3D();
	virt_elastic_matrix3d[2] = new ElasticMatrix3D();
	random_axis = NULL;
	random_axis_inv = NULL;
	basis_quantity = 0;
	volume_calc = new SimpleVolumeCalculator();
	free_border_calc = new FreeBorderCalculator();
	fixed_border_calc = new FixedBorderCalculator();
	ext_force_calc = new ExternalForceCalculator();
	ext_v_calc = new ExternalVelocityCalculator();
	ext_val_calc = new ExternalValuesCalculator();
	adhesion_contact_calc = new AdhesionContactCalculator();
	sph_connector = new SPHConnector();
};

GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::~GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis()
{
	delete(elastic_matrix3d[0]);
	delete(elastic_matrix3d[1]);
	delete(elastic_matrix3d[2]);
	delete(virt_elastic_matrix3d[0]);
	delete(virt_elastic_matrix3d[1]);
	delete(virt_elastic_matrix3d[2]);
	if(random_axis != NULL)
		free(random_axis);
	if(random_axis_inv != NULL)
		free(random_axis_inv);
};

void GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::prepare_part_step(ElasticNode* cur_node, ElasticMatrix3D* matrix, int stage, int basis_num)
{

	if(stage < 3) {
		matrix->prepare_matrix( cur_node->la, cur_node->mu, cur_node->rho, random_axis_inv[basis_num].ksi[0][stage], 
				random_axis_inv[basis_num].ksi[1][stage], random_axis_inv[basis_num].ksi[2][stage] );
	} else {
		throw GCMException( GCMException::METHOD_EXCEPTION, "Bad stage number");
	}
};

void GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::drop_deviator(ElasticNode* cur_node, ElasticNode* new_node)
{

	float J2 = sqrt( ( (cur_node->values[3] - cur_node->values[6]) * (cur_node->values[3] - cur_node->values[6]) 
			+ (cur_node->values[6] - cur_node->values[8]) * (cur_node->values[6] - cur_node->values[8])
			+ (cur_node->values[3] - cur_node->values[8]) * (cur_node->values[3] - cur_node->values[8])
			+ 6 * ( (cur_node->values[4]) * (cur_node->values[4]) + (cur_node->values[5]) * (cur_node->values[5])
				+ (cur_node->values[7]) * (cur_node->values[7])) ) / 6 );

	float p = (cur_node->values[3] + cur_node->values[6] + cur_node->values[8]) / 3;

	if (cur_node->yield_limit < J2)
	{
		new_node->values[3] = ((cur_node->values[3] - p) * cur_node->yield_limit / J2) + p;
		new_node->values[4] = cur_node->values[4] * cur_node->yield_limit / J2;
		new_node->values[5] = cur_node->values[5] * cur_node->yield_limit / J2;
		new_node->values[6] = ((cur_node->values[6] - p) * cur_node->yield_limit / J2) + p;
		new_node->values[7] = cur_node->values[7] * cur_node->yield_limit / J2;
		new_node->values[8] = ((cur_node->values[8] - p) * cur_node->yield_limit / J2) + p;
	}
};

int GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::prepare_node(ElasticNode* cur_node, ElasticMatrix3D* matrixes[], float time_step, int stage, TetrMesh* mesh, float dksi[], bool inner[], ElasticNode previous_nodes[], float outer_normal[], int ppoint_num[], int basis_num)
{
	return prepare_node(cur_node, matrixes, time_step, stage, mesh, dksi, inner, previous_nodes, outer_normal, ppoint_num, basis_num, false);
};

int GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::prepare_node(ElasticNode* cur_node, ElasticMatrix3D* matrixes[], float time_step, int stage, TetrMesh* mesh, float dksi[], bool inner[], ElasticNode previous_nodes[], float outer_normal[], int ppoint_num[], int basis_num, bool debug)
{
	if(debug)
		*logger < "DEBUG 1";

	if (cur_node->isBorder ())
		mesh->find_border_node_normal(cur_node->local_num, &outer_normal[0], &outer_normal[1], &outer_normal[2]);

	if(debug)
		*logger < "DEBUG 2";

	//  Prepare matrixes  A, Lambda, Omega, Omega^(-1)
	prepare_part_step(cur_node, matrixes[stage], stage, basis_num);

	if(debug)
		*logger < "DEBUG 3";

	for(int i = 0; i < 9; i++)
		dksi[i] = - matrixes[stage]->L(i,i) * time_step;

	float alpha = 0;

	int outer_count;

	if(debug)
		*logger < "DEBUG 4";

	do {
		outer_count = find_nodes_on_previous_time_layer(cur_node, stage, mesh, alpha, dksi, inner, previous_nodes, outer_normal, ppoint_num, basis_num, debug);
		alpha += 0.1;
	} while( (outer_count != 0) && (outer_count != 3) && (alpha < 1.01) );

	return outer_count;
};

void GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::log_node_diagnostics(ElasticNode* cur_node, int stage, float outer_normal[], TetrMesh* mesh, int basis_num, ElasticMatrix3D* matrixes[], float time_step, ElasticNode previous_nodes[], int ppoint_num[], bool inner[], float dksi[])
{
	int outer_count = prepare_node(cur_node, matrixes, time_step, stage, mesh, dksi, inner, previous_nodes, outer_normal, ppoint_num, basis_num, true);

	*logger < "NODE DETAILS:";
	*logger << "STAGE: " < stage;
	*logger << "NODE " << cur_node->local_num << ": x: " << cur_node->coords[0] 
				<< " y: " << cur_node->coords[1]
				<< " z: " < cur_node->coords[2];
	*logger << "OUTER COUNT: " < outer_count;
	*logger < "VALUES: ";
	for(int j = 0; j < 9; j++)
		*logger << "Value[" << j << "] = " < cur_node->values[j];
	if( cur_node->isBorder ()) {
		*logger < "BORDER";
		if( ( cur_node->contact_data != NULL ) && ( cur_node->contact_data->axis_plus[stage] == -1 ) && ( cur_node->contact_data->axis_minus[stage] == -1 ) ) {
			*logger < "BORDER WITHOUT CONTACT";
		} else {
			*logger < "CONTACT BORDER";

			if( cur_node->contact_data != NULL ) {
				*logger < "CONTACT DATA:";
				for(int k = 0; k < 3; k++)
					*logger << "Axis[" << k << "]: minus: " << cur_node->contact_data->axis_minus[k] << " plus: " < cur_node->contact_data->axis_plus[k];

				ElasticNode* virt_node;
				if( cur_node->contact_data->axis_plus[stage] != -1 )
					virt_node=mesh->mesh_set->getNode(cur_node->contact_data->axis_plus[stage]);
				else
					virt_node=mesh->mesh_set->getNode(cur_node->contact_data->axis_minus[stage]);

				*logger << "VIRT NODE " << virt_node->local_num << ":"
						<< " x: " << virt_node->coords[0]
						<< " y: " << virt_node->coords[1]
						<< " z: " < virt_node->coords[2];
				*logger < "VIRT NODE VALUES: ";
				for(int j = 0; j < 9; j++)
					*logger << "Value[" << j << "] = " < virt_node->values[j];

			}
		}
	} else {
		*logger < "INNER";
	}
	*logger << "OUTER_NORMAL: " << outer_normal[0] << " " << outer_normal[1] << " " < outer_normal[2];
	*logger << "NEIGH: " < (cur_node->elements)->size();
	for(int i = 0; i < (cur_node->elements)->size(); i++) {
		Tetrahedron* tmp_tetr = mesh->get_tetrahedron( (cur_node->elements)->at(i) );
		*logger << "\tTetr " << tmp_tetr->local_num << " Neigh_num: " < i;
		for(int j = 0; j < 4; j++) {
			ElasticNode* tmp_node = mesh->get_node( tmp_tetr->vert[j] );
			*logger << "\t\tVert: " << j << " num: " << tmp_node->local_num << "\t"
					<< " x: " << tmp_node->coords[0]
					<< " y: " << tmp_node->coords[1]
					<< " z: " < tmp_node->coords[2];
		}
	}
	for(int i = 0; i < 3; i++)
		*logger << "KSI[" << i << "]: x: " << random_axis[basis_num].ksi[i][0]
				<< " y: " << random_axis[basis_num].ksi[i][1]
				<< " z: " < random_axis[basis_num].ksi[i][2];
	for(int i = 0; i < 9; i++) {
		if(inner[i]) {
			*logger << "INNER OMEGA: num: " << i
				<< " val: " << matrixes[stage]->L(i,i)
				<< " step: " < matrixes[stage]->L(i,i) * time_step;
		} else {
			*logger << "OUTER OMEGA: num: " << i 
				<< " val: " << matrixes[stage]->L(i,i)
				<< " step: " < matrixes[stage]->L(i,i) * time_step;
		}
		*logger << "\t Point x: " << previous_nodes[ppoint_num[i]].coords[0]
				<< " y: " << previous_nodes[ppoint_num[i]].coords[1]
				<< " z: " < previous_nodes[ppoint_num[i]].coords[2];
	}
};

void GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::do_next_part_step(ElasticNode* cur_node, ElasticNode* new_node, float time_step, int stage, TetrMesh* mesh)
{
	for(int i = 0; i < 9; i++)
		if( isnan(cur_node->values[i]) || isinf(cur_node->values[i]) )
			throw GCMException(GCMException::METHOD_EXCEPTION, "NAN or INF detected. Instability?");

	if(stage == 3) {
		drop_deviator(cur_node, new_node);
		return;
	}

	// Variables used in calculations internally

	// Delta x on previous time layer for all the omegas
	// 	omega_new_time_layer(ksi) = omega_old_time_layer(ksi+dksi)
	float dksi[9];

	// If the corresponding point on previous time layer is inner or not
	bool inner[9];

	// We will store interpolated nodes on previous time layer here
	// We know that we need five nodes for each direction (corresponding to Lambdas -C1, -C2, 0, C2, C1)
	// TODO  - We can  deal with (lambda == 0) separately
	ElasticNode previous_nodes[5];

	// Outer normal at current point
	float outer_normal[3];

	// This array will link omegas with corresponding interpolated nodes they should be copied from
	int ppoint_num[9];

	int basis_num = cur_node->local_num;

	// Here we will store (omega = Matrix_OMEGA * u)
	float omega[9];

	// Number of outer characteristics
	int outer_count = prepare_node(cur_node, elastic_matrix3d, time_step, stage, mesh, dksi, inner, previous_nodes, outer_normal, ppoint_num, basis_num);

	float* previous_values[9];
	for(int i = 0; i < 9; i++)
		previous_values[i] = &previous_nodes[ppoint_num[i]].values[0];

	// TODO - merge this condition with the next ones
	if((outer_count != 0) && (outer_count != 3)) {
		*logger << "There are " << outer_count < " 'outer' characteristics for real node.";
		outer_count = prepare_node(cur_node, elastic_matrix3d, time_step, stage, mesh, dksi, inner, previous_nodes, outer_normal, ppoint_num, basis_num, true);
		log_node_diagnostics(cur_node, stage, outer_normal, mesh, basis_num, elastic_matrix3d, time_step, previous_nodes, ppoint_num, inner, dksi);
		throw GCMException(GCMException::METHOD_EXCEPTION, "Illegal number of outer characteristics");
	}

	// If all the omegas are 'inner'
	// omega = Matrix_OMEGA * u
	// new_u = Matrix_OMEGA^(-1) * omega
	// TODO - to think - if all omegas are 'inner' can we skip matrix calculations and just use new_u = interpolated_u ?
	if( outer_count == 0 )
	{
		// Special cases - smth bad happens
		// TODO - commenting it out - it is possible that border node is inner for some axis
		// if( cur_node->border_type == BORDER )	// node is marked as border
		// {
		// 	if(logger != NULL)
		// 	log_node_diagnostics(cur_node, stage, outer_normal, mesh, basis_num, elastic_matrix3d, time_step, previous_nodes, ppoint_num, inner, dksi);
		// 	throw GCMException(GCMException::METHOD_EXCEPTION, "outer_count == 0 for BORDER node");
		// }

		cur_node->volume_calculator->do_calc(new_node, elastic_matrix3d[stage], previous_values);

	}
	// If there are 3 'outer' omegas - we should use border or contact algorithm
	// TODO - ... we should also use it when l*t/h > 1 and there is 1 'outer' omega
	// 		(we should add 2 more corresponding omegas to 'outer' manually
	else if ( outer_count == 3 )
	{
		// Check contact state

		// Special cases - smth bad happens
		// Node is not border at all
		if ( !cur_node->isBorder ())
		{
			log_node_diagnostics(cur_node, stage, outer_normal, mesh, basis_num, elastic_matrix3d, time_step, previous_nodes, ppoint_num, inner, dksi);
			throw GCMException(GCMException::METHOD_EXCEPTION, "Border node is not marked as border");
		}
		// Node should be border but it has no contact_data struct
		else if( cur_node->contact_data == NULL )
		{
			log_node_diagnostics(cur_node, stage, outer_normal, mesh, basis_num, elastic_matrix3d, time_step, previous_nodes, ppoint_num, inner, dksi);
			throw GCMException(GCMException::METHOD_EXCEPTION, "Border node has no contact data struct");
		}
		// Both directions are marked as contact
		else if ( ( cur_node->contact_data->axis_plus[stage] != -1 ) && ( cur_node->contact_data->axis_minus[stage] != -1 ) )
		{
			log_node_diagnostics(cur_node, stage, outer_normal, mesh, basis_num, elastic_matrix3d, time_step, previous_nodes, ppoint_num, inner, dksi);
			throw GCMException(GCMException::METHOD_EXCEPTION, "Both direction of contact are marked as contact");
		}

		// If both directions show no contact - use border algorithm, otherwise use contact algorithm
		if( ( cur_node->contact_data->axis_plus[stage] == -1 ) && ( cur_node->contact_data->axis_minus[stage] == -1 ) )
		{
			// We use basis axis here instead of outer_normal because:
			//    - for 'normal' points first axis coincides with normal and that's it
			//    - for edges and verts it is the only way to avoid singular matrix
			// Singular matrix appears because:
			//    - current axis is used for A calculation, thus for 6 equations in Omega
			//    - outer normal gives us 3 additional equations to replace Omega's ones
			//    - normal has no projection on axis for all axis except the first.
			// Effectively this approach 'smooth' edges and verts.
			cur_node->border_condition->do_calc(mesh->get_current_time(), cur_node->coords, new_node, elastic_matrix3d[stage], previous_values, inner, (random_axis + basis_num)->ksi[stage]);

		} else {

			ElasticNode* virt_node;
			if( cur_node->contact_data->axis_plus[stage] != -1 )
				virt_node = mesh->mesh_set->getNode( cur_node->contact_data->axis_plus[stage] );
			else 
				virt_node = mesh->mesh_set->getNode( cur_node->contact_data->axis_minus[stage] );

			// Mark virt node as having contact state
			// TODO FIXME - most probably CollisionDetector should do it
			// But we should check it anycase
			virt_node->contact_data = (contact_state*) malloc(sizeof(contact_state));
			mesh->clear_contact_data(virt_node);
			if( cur_node->contact_data->axis_plus[stage] != -1 )
				virt_node->contact_data->axis_minus[stage] = cur_node->contact_data->axis_plus[stage];
			else
				virt_node->contact_data->axis_plus[stage] = cur_node->contact_data->axis_minus[stage];

			// Variables used in calculations internally

			// Delta x on previous time layer for all the omegas
			// 	omega_new_time_layer(ksi) = omega_old_time_layer(ksi+dksi)
			float virt_dksi[9];

			// If the corresponding point on previous time layer is inner or not
			bool virt_inner[9];

			// We will store interpolated nodes on previous time layer here
			// We know that we need five nodes for each direction (corresponding to Lambdas -C1, -C2, 0, C2, C1)
			// TODO  - We can  deal with (lambda == 0) separately
			ElasticNode virt_previous_nodes[5];

			// Outer normal at current point
			float virt_outer_normal[3];

			// This array will link omegas with corresponding interpolated nodes they should be copied from
			int virt_ppoint_num[9];

			// Here we will store (omega = Matrix_OMEGA * u)
			float virt_omega[9];

			// Number of outer characteristics
			int virt_outer_count = prepare_node(virt_node, virt_elastic_matrix3d, time_step, stage, virt_node->mesh, virt_dksi, virt_inner, virt_previous_nodes, virt_outer_normal, virt_ppoint_num, basis_num);

			float* virt_previous_values[9];
			for(int i = 0; i < 9; i++)
				virt_previous_values[i] = &virt_previous_nodes[virt_ppoint_num[i]].values[0];

			// TODO - merge this condition with the next ones
			if( virt_outer_count != 3 ) {
					*logger << "There are " << virt_outer_count < " 'outer' characteristics for virt node.";
					*logger << "MESH " << mesh->zone_num << " REAL NODE " << cur_node->local_num << ": " 
								<< "x: " << cur_node->coords[0] 
								<< " y: " << cur_node->coords[1] 
								<< " z: " < cur_node->coords[2];
				log_node_diagnostics(virt_node, stage, virt_outer_normal, virt_node->mesh, basis_num, virt_elastic_matrix3d, time_step, virt_previous_nodes, virt_ppoint_num, virt_inner, virt_dksi);
				throw GCMException(GCMException::METHOD_EXCEPTION, "Illegal number of outer characteristics");
			}

			// Check that 'paired node' is in the direction of 'outer' characteristics
			// If it is not the case - we have strange situation when 
			// we replace 'outer' points data with data of 'paired node' from different axis direction.
			
			// For all characteristics of real node and virt node
			for(int i = 0; i < 9; i++)
			{
				float v_x_outer[3];
				float v_x_virt[3];
				// Real node - if characteristic is 'outer'
/*				if(!inner[i])
				{
					// Find directions to corresponding 'outer' point and to virt 'paired node'
					for(int j = 0; j < 3; j++) {
						v_x_outer[j] = previous_nodes[ppoint_num[i]].coords[j] - cur_node->coords[j];
						v_x_virt[j] = virt_node->coords[j] - cur_node->coords[j];
					}
					// If directions are different - smth bad happens
					if( (v_x_outer[0] * v_x_virt[0]
						 + v_x_outer[1] * v_x_virt[1] + v_x_outer[2] * v_x_virt[2]) < 0 )
					{
						*logger << "MESH " << mesh->zone_num << "REAL NODE " << cur_node->local_num << ": " 
								<< "x: " << cur_node->coords[0] 
								<< " y: " << cur_node->coords[1] 
								<< " z: " < cur_node->coords[2];
						log_node_diagnostics(cur_node, stage, outer_normal, mesh, basis_num, elastic_matrix3d, time_step, previous_nodes, ppoint_num, inner, dksi);
						*logger << "'Outer' direction: " << v_x_outer[0] << " " 
							<< v_x_outer[1] << " " < v_x_outer[2];
						*logger << "'Virt' direction: " << v_x_virt[0] << " "
							<< v_x_virt[1] << " " < v_x_virt[2];
						throw GCMException( GCMException::METHOD_EXCEPTION, "Bad contact from real node point of view: 'outer' and 'virt' directions are different");
					}
				}*/
// We switch it off because it conflicts sometimes with 'safe_direction'
/*				// Virt node - if characteristic is 'outer'
				if(!virt_inner[i])
				{
					// Find directions to corresponding 'outer' point and to real 'paired node'
					for(int j = 0; j < 3; j++) {
						v_x_outer[j] = virt_previous_nodes[virt_ppoint_num[i]].coords[j] - virt_node->coords[j];
						v_x_virt[j] = cur_node->coords[j] - virt_node->coords[j];
					}
					// If directions are different - smth bad happens
					if( (v_x_outer[0] * v_x_virt[0]
						+ v_x_outer[1] * v_x_virt[1] + v_x_outer[2] * v_x_virt[2]) < 0 )
					{
						*logger << "MESH " << mesh->zone_num << "REAL NODE " << cur_node->local_num << ": " 
								<< "x: " << cur_node->coords[0] 
								<< " y: " << cur_node->coords[1] 
								<< " z: " < cur_node->coords[2];
						log_node_diagnostics(virt_node, stage, virt_outer_normal, virt_node->mesh, basis_num, virt_elastic_matrix3d, time_step, virt_previous_nodes, virt_ppoint_num, virt_inner, virt_dksi);
						*logger << "'Outer' direction: " << v_x_outer[0] << " "
							<< v_x_outer[1] << " "< v_x_outer[2];
						*logger << "'Virt' direction: " << v_x_virt[0] << " "
							<< v_x_virt[1] << " " < v_x_virt[2];
						throw GCMException( GCMException::METHOD_EXCEPTION, "Bad contact from virt node point of view: 'outer' and 'virt' directions are different");
					}
				}*/
			}

			cur_node->contact_condition->do_calc(mesh->get_current_time(), cur_node, virt_node, new_node, elastic_matrix3d[stage], previous_values, inner, virt_elastic_matrix3d[stage], virt_previous_values, virt_inner, outer_normal);

			free(virt_node->contact_data);
		}
	}
};

int GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::find_nodes_on_previous_time_layer(ElasticNode* cur_node, int stage, TetrMesh* mesh, float alpha, float dksi[], bool inner[], ElasticNode previous_nodes[], float outer_normal[], int ppoint_num[], int basis_num)
{
	return find_nodes_on_previous_time_layer(cur_node, stage, mesh, alpha, dksi, inner, previous_nodes, outer_normal, ppoint_num, basis_num, false);
};

int GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::find_nodes_on_previous_time_layer(ElasticNode* cur_node, int stage, TetrMesh* mesh, float alpha, float dksi[], bool inner[], ElasticNode previous_nodes[], float outer_normal[], int ppoint_num[], int basis_num, bool debug)
{
	float safe_direction[3];
	// FIXME - bad choice - just for PoC
	for(int j = 0; j < 3; j++)
		safe_direction[j] = - ( cur_node->coords[j] - (mesh->nodes).at(cur_node->local_num).coords[j] );

	int tetr_num = cur_node->elements->at(0);
	Tetrahedron* neighTetr = mesh->get_tetrahedron(tetr_num);

	float tetr_center[3];
	// Tetr center
	for(int i = 0; i < 3; i++)
		tetr_center[i] = ( (mesh->nodes)[ neighTetr->vert[0] ].coords[i] + (mesh->nodes)[ neighTetr->vert[1] ].coords[i]
				+ (mesh->nodes)[ neighTetr->vert[2] ].coords[i] + (mesh->nodes)[ neighTetr->vert[3] ].coords[i] ) / 4;

//	for(int i = 0; i < 3; i++)
//		safe_direction[i] = - ( cur_node->coords[i] - tetr_center[i] );


	if(debug)
		*logger << "DEBUG 5 alpha=" < alpha;


	if( (alpha > 1.01) || (alpha < 0) )
		throw GCMException( GCMException::METHOD_EXCEPTION, "Bad alpha");

	if (stage >= 3)
		throw GCMException( GCMException::METHOD_EXCEPTION, "Bad stage number");

	if(alpha > 1.0) {
		alpha = 1.0;
//		*logger < "Shitty alpha";
	}

	// Just tmp tetr pointer
	Tetrahedron* tmp_tetr;

	int count = 0;

	float dx[3];
	float dx_ksi[3];
	float safe_direction_projection[3];

	// For all omegas
	for(int i = 0; i < 9; i++)
	{
		// Check prevoius omegas ...
		bool already_found = false;
		for(int j = 0; j < i; j++)
		{
			// ... And try to find if we have already worked with the required point
			// on previous time layer (or at least with the point that is close enough)
			if( fabs(dksi[i] - dksi[j]) <= 0.01 * fabs(dksi[i] + dksi[j]) ) // TODO - avoid magick number!
			{
				// If we have already worked with this point - just remember the number
				already_found = true;
				ppoint_num[i] = ppoint_num[j];
				inner[i] = inner[j];
			}
		}
		// If we do not have necessary point in place - ...
		if(already_found == false)
		{
			// ... Put new number ...
			ppoint_num[i] = count;

			previous_nodes[count] = *cur_node;

			// ... Find vectors ...
			for(int j = 0; j < 3; j++)
				dx_ksi[j] = dksi[i] * random_axis[basis_num].ksi[stage][j];

			// TODO - this criteria is used when we use rotation_matrix_around_normal
			// In this case dx_ksi_normal_projection_modul != 0
			// float dx_ksi_normal_projection_modul = dx_ksi[0] * outer_normal[0] + dx_ksi[1] * outer_normal[1] + dx_ksi[2] * outer_normal[2];
			// TODO - this criteria to be used with rotation_matrix_with_normal
			// float dx_ksi_normal_projection_modul = - fabs(dksi[i]);

			// FIXME - bad choice - just for PoC
			float safe_direction_projection_modul = 1;//0.95 * ((TetrMesh_1stOrder*)mesh)->get_min_h();

			// ... Calculate coordinates ...
			for(int j = 0; j < 3; j++) {
				// If we have contact virt 'paired node' -
				//  we should alter direction without contact and guarantee that
				//  exactly this direction will give us internal node after all
				if( ( cur_node->contact_data != NULL ) 
					&& ( ( cur_node->contact_data->axis_plus[stage] != -1 ) 
						|| ( cur_node->contact_data->axis_minus[stage] != -1 ) ) ) {
					// Positive direction has virtual node - alter negative direction only
					// and use opposite direction to outer normal as last chance
					if( cur_node->contact_data->axis_plus[stage] != -1 )
					{
						if( dksi[i] >= 0 )
							dx[j] = dx_ksi[j];
						else
						{
							safe_direction_projection[j] = safe_direction_projection_modul * safe_direction[j];
							dx[j] = dx_ksi[j] * (1 - alpha) + safe_direction_projection[j] * alpha;
						}
					}
					// Negative direction has virtual node - alter positive direction only
					// and use opposite direction to outer normal as last chance
					else if ( cur_node->contact_data->axis_minus[stage] != -1 )
					{
						if( dksi[i] <= 0 )
							dx[j] = dx_ksi[j];
						else
						{
							safe_direction_projection[j] = safe_direction_projection_modul * safe_direction[j];
							dx[j] = dx_ksi[j] * (1 - alpha) + safe_direction_projection[j] * alpha;
						}
					}
					else
					{
						throw GCMException( GCMException::METHOD_EXCEPTION, "Bad contact data");
					}
				// Otherwise (internal node or border node without contact) - 
				//  alter both directions and wait for one of them to give internal node
				} else {
					// FIXME - we alter only one direction (chosen at will for the moment) 
					// to avoid case when both can become internal and give us outer_count == 0 (!)
					if( dksi[i] >= 0 ) {
						dx[j] = dx_ksi[j];
					} else {
						safe_direction_projection[j] = safe_direction_projection_modul * safe_direction[j];
						dx[j] = dx_ksi[j] * (1 - alpha) + safe_direction_projection[j] * alpha;
					}
				}

				// This difference is not zero when cur_node coords were altered -
				//  	it happens if it is virtual node moved compared with its former base one
				dx[j] += ( cur_node->coords[j] - (mesh->nodes).at(cur_node->local_num).coords[j] );

				previous_nodes[count].coords[j] = (mesh->nodes).at(cur_node->local_num).coords[j] + dx[j];
			}

			// FIXME last chance attempt
			if(alpha > 0.95) {
				float fafa = ((TetrMesh_1stOrder*)mesh)->get_min_h();
//*logger < "Last attempt happens";
//for(int z = 0; z < 3; z++)
	//*logger < cur_node->coords[z];
//				if( dksi[i] < 0 )
				if( ! ( (cur_node->isBorder ()) && (stage == 0) && (dksi[i] > 0) ) )
					for(int z = 0; z < 3; z++) {
//						dx[z] = - fafa * safe_direction[z];
						dx[z] = tetr_center[z] - (mesh->nodes).at(cur_node->local_num).coords[z];
						previous_nodes[count].coords[z] = (mesh->nodes).at(cur_node->local_num).coords[z] + dx[z];
//						*logger < dx[z];
					}
//			dx[z] = safe_direction[z] + ( cur_node->coords[z] - (mesh->nodes).at(cur_node->local_num).coords[z] );
//			dx[z] = tetr_center[z] - (mesh->nodes).at(cur_node->local_num).coords[z];
			}

			if(debug)
				*logger < "DEBUG 6";

			// ... Find owner tetrahedron ...
			tmp_tetr = mesh->find_owner_tetr(cur_node, dx[0], dx[1], dx[2], debug);

			// special cases - contact direction must be outer - we have issues with virt nodes when they are moved
			// and some characterictics (shorter ones) DO fit into inner area
			if ( ( cur_node->contact_data != NULL ) 
					&& ( cur_node->contact_data->axis_plus[stage] != -1 ) && ( dksi[i] > 0 ) )
			{
				inner[i] = false;
			}
			else if ( ( cur_node->contact_data != NULL ) 
					&& ( cur_node->contact_data->axis_minus[stage] != -1 ) && ( dksi[i] < 0 ) )
			{
				inner[i] = false;
			}
			else if( tmp_tetr != NULL )
			{
				if(debug)
					*logger < "DEBUG 6.1";
				// ... And interpolate values
				mesh->interpolate(&previous_nodes[count], tmp_tetr);
				inner[i] = true;
			}
			else if (dksi[i] == 0)
			{
				if(debug)
					*logger < "DEBUG 6.2";
				previous_nodes[count] = *cur_node;
				inner[i] = true;
			}
			else if ( (cur_node->isInner ()) /*|| (stage != 0)*/ )
			{
				*logger < "We need new method here!";
				*logger << cur_node->local_num << " " << cur_node->coords[0] << " "
						 << cur_node->coords[1] << " " < cur_node->coords[2];
				for(int j = 0; j < 3; j++)
					*logger < dksi[i] * random_axis[basis_num].ksi[stage][j];

				ElasticNode cross;
				tmp_tetr = mesh->find_border_cross(cur_node, dx[0], dx[1], dx[2], &cross);

				for(int j = 0; j < 3; j++)
					previous_nodes[count].coords[j] = cross.coords[j];

				if( ! cross.isOwnedBy( SPH ) )
				{
					// !!!! FIXME - time-aware interpolation required
					mesh->interpolate(&previous_nodes[count], tmp_tetr);
				}
				else
				{
					sph_connector->interpolate( previous_nodes[count].coords, previous_nodes[count].values );
				}
				inner[i] = true;
			}
			else if ( cur_node->isOwnedBy(SPH) )
			{
				sph_connector->interpolate( previous_nodes[count].coords, previous_nodes[count].values );
				inner[i] = true;
			}
			else
			{
				if(debug) {
					*logger < "DEBUG 6.3";
					float fafa = ((TetrMesh_1stOrder*)mesh)->get_min_h();
					if( mesh->find_owner_tetr(cur_node, 
						-fafa * safe_direction[0], -fafa * safe_direction[1], -fafa * safe_direction[2], debug) != NULL )
					{
						*logger << "WTF? safe_direction is correct. Modul is " 
								<< fabs(safe_direction_projection_modul) << " Safe value is: " < fafa;
						for(int j = 0; j < 3; j++)
							*logger << "DATA " << dx[j] << " " << dx_ksi[j] << " " 
									<< safe_direction_projection[j] << " " 
									<< (cur_node->coords[j] - (mesh->nodes).at(cur_node->local_num).coords[j]) 
									<< " " < alpha;
					} else
						*logger < "Normal is really incorrect";
				}
				// There is no value in interpolation in this case
				//	as long as characteristic is out of region
				inner[i] = false;
			}
			count++;
		}
	}

	int outer_count = 0;
	for(int i = 0; i < 9; i++)
		if(!inner[i])
			outer_count++;

	return outer_count;
};

int GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::get_number_of_stages()
{
	return 4;
};

float GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::get_max_lambda(ElasticNode* node, TetrMesh* mesh)
{
	create_random_axis(node, mesh);

	// We just return sqrt((la+2*mu)/rho) because axis are randomized by rotation, so x^2+y^2+z^2 == 1
	// TODO - explicit check?
	return sqrt( ( (node->la) + 2 * (node->mu) ) / (node->rho) );
};

int GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::create_random_axis(ElasticNode* cur_node, TetrMesh* mesh)
{

	if(basis_quantity < (mesh->nodes).size()) {
		basis_quantity = (mesh->nodes).size();
		random_axis = (basis*)realloc(random_axis, sizeof(basis) * basis_quantity);
		random_axis_inv = (basis*)realloc(random_axis_inv, sizeof(basis) * basis_quantity);
	}

	int basis_num = cur_node->local_num;
	const double PI = atan(1) * 4;

	// Outer normal at current point
	float outer_normal[3];

	if(cur_node->isInner ()) {

		// TODO think about limits - PI or 2*PI
		float phi = (rand() % 360) * 2 * PI / 360;
		float psi = (rand() % 360) * 2 * PI / 360;
		float teta = (rand() % 360) * 2 * PI / 360;

		create_rotation_matrix(cur_node->local_num, phi, psi, teta);

		// FIXME
		create_E_matrix(cur_node->local_num);

	} else if (cur_node->isBorder ()) {

		mesh->find_border_node_normal(cur_node->local_num, &outer_normal[0], &outer_normal[1], &outer_normal[2]);

		// TODO think about limits - PI or 2*PI
		float teta = (rand() % 360) * 2 * PI / 360;

		// Rotate at random angle around normal
		create_rotation_matrix_with_normal(cur_node->local_num, outer_normal[0], outer_normal[1], outer_normal[2], teta);

	} else {
		throw GCMException( GCMException::METHOD_EXCEPTION, "Unknown border type");
	}

	// FIXME
	// create_E_matrix(cur_node->local_num);

	// Attach new basis to node - we need it for CollisionDetector to create virtual nodes using directions of these axis
	// TODO - should it be this way at all?
	cur_node->local_basis = &random_axis[basis_num];

	find_inverse_matrix(basis_num);

	return 0;
};

// FIXME
void GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::create_E_matrix(int node_num)
{
	random_axis[node_num].ksi[0][0] = 1;
	random_axis[node_num].ksi[0][1] = 0;
	random_axis[node_num].ksi[0][2] = 0;

	random_axis[node_num].ksi[1][0] = 0;
        random_axis[node_num].ksi[1][1] = 1;
        random_axis[node_num].ksi[1][2] = 0;

	random_axis[node_num].ksi[2][0] = 0;
        random_axis[node_num].ksi[2][1] = 0;
        random_axis[node_num].ksi[2][2] = 1;
};

void GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::create_rotation_matrix(int node_num, float phi, float psi, float teta)
{
	random_axis[node_num].ksi[0][0] = cos(teta) * cos(psi);
	random_axis[node_num].ksi[0][1] = cos(teta) * sin(psi);
	random_axis[node_num].ksi[0][2] = - sin(teta);

	random_axis[node_num].ksi[1][0] = - cos(phi) * sin(psi) + sin(phi) * sin(teta) * cos(psi);
	random_axis[node_num].ksi[1][1] = cos(phi) * cos(psi) + sin(phi) * sin(teta) * sin(psi);
	random_axis[node_num].ksi[1][2] = sin(phi) * cos(teta);

	random_axis[node_num].ksi[2][0] = sin(phi) * sin(psi) + cos(phi) * sin(teta) * cos(psi);
	random_axis[node_num].ksi[2][1] = - sin(phi) * cos(psi) + cos(phi) * sin(teta) * sin(psi);
	random_axis[node_num].ksi[2][2] = cos(phi) * cos(teta);
};

void GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::create_rotation_matrix_around_normal(int node_num, float x, float y, float z, float teta)
{
	float sqrt3 = sqrt(3);
	float sqrt2 = sqrt(2);
	float sqrt6 = sqrt2 * sqrt3;

	basis new_axis;

	// Rotation from normal "quasi-(1;1;1)" to first axis "quasi-(1;0;0)"
	float x_m[3][3];
	x_m[0][0] = sqrt3/3;	x_m[1][0] = sqrt3/3;	x_m[2][0] = sqrt3/3;
	x_m[0][1] = -sqrt2/2;	x_m[1][1] = sqrt2/2;	x_m[2][1] = 0;
	x_m[0][2] = -sqrt6/6;	x_m[1][2] = -sqrt6/6;	x_m[2][2] = 2*sqrt6/6;

	for(int i = 0; i < 3; i++)
		new_axis.ksi[0][i] = x_m[0][i] * x + x_m[1][i] * y + x_m[2][i] * z;

	// Rotation from normal "quasi-(1;1;1)" to second axis "quasi-(0;1;0)"
	float y_m[3][3];
	y_m[0][0] = sqrt2/2;	y_m[1][0] = -sqrt2/2;	y_m[2][0] = 0;
	y_m[0][1] = sqrt3/3;	y_m[1][1] = sqrt3/3;	y_m[2][1] = sqrt3/3;
	y_m[0][2] = -sqrt6/6;	y_m[1][2] = -sqrt6/6;	y_m[2][2] = 2*sqrt6/6;

	for(int i = 0; i < 3; i++)
		new_axis.ksi[1][i] = y_m[0][i] * x + y_m[1][i] * y + y_m[2][i] * z;

	// Rotation from normal "quasi-(1;1;1)" to third axis "quasi-(0;0;1)"
	float z_m[3][3];
	z_m[0][0] = sqrt2/2;	z_m[1][0] = 0;		z_m[2][0] = -sqrt2/2;
	z_m[0][1] = -sqrt6/6;	z_m[1][1] = 2*sqrt6/6;	z_m[2][1] = -sqrt6/6;
	z_m[0][2] = sqrt3/3;	z_m[1][2] = sqrt3/3;	z_m[2][2] = sqrt3/3;

	for(int i = 0; i < 3; i++)
		new_axis.ksi[2][i] = z_m[0][i] * x + z_m[1][i] * y + z_m[2][i] * z;

	// Rotation matrix - random teta around normal
	float rot[3][3];

	rot[0][0] = cos(teta) + x * x * (1 - cos(teta));
	rot[0][1] = y * x * (1 - cos(teta)) + z * sin(teta);
	rot[0][2] = z * x * (1 - cos(teta)) - y * sin(teta);

	rot[1][0] = y * x * (1 - cos(teta)) - z * sin(teta);
	rot[1][1] = cos(teta) + y * y * (1 - cos(teta));
	rot[1][2] = y * z * (1 - cos(teta)) + x * sin(teta);

	rot[2][0] = x * z * (1 - cos(teta)) + y * sin(teta);
	rot[2][1] = y * z * (1 - cos(teta)) - x * sin(teta);
	rot[2][2] = cos(teta) + z * z * (1 - cos(teta));

	// Rotate
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			random_axis[node_num].ksi[i][j] = 0;
			for(int k = 0; k < 3; k++) {
				random_axis[node_num].ksi[i][j] += rot[k][j] * new_axis.ksi[i][k];
			}
		}
	}

};

void GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::create_rotation_matrix_with_normal(int node_num, float x, float y, float z, float teta)
{
	basis new_axis;

	// First axis goes along normal
	new_axis.ksi[0][0] = x;
	new_axis.ksi[0][1] = y;
	new_axis.ksi[0][2] = z;

	// Second axis is perpendicular to it
	if( (fabs(x) >= fabs(y)) && (fabs(x) >= fabs(z)) ) {
		new_axis.ksi[1][0] = y;
		new_axis.ksi[1][1] = -x;
		new_axis.ksi[1][2] = 0;
	} else if( (fabs(y) >= fabs(x)) && (fabs(y) >= fabs(z)) ) {
		new_axis.ksi[1][0] = y;
		new_axis.ksi[1][1] = -x;
		new_axis.ksi[1][2] = 0;
	} else {
		new_axis.ksi[1][0] = 0;
		new_axis.ksi[1][1] = -z;
		new_axis.ksi[1][2] = y;
	}
	float modul = sqrt( new_axis.ksi[1][0] * new_axis.ksi[1][0] + new_axis.ksi[1][1] * new_axis.ksi[1][1] + new_axis.ksi[1][2] * new_axis.ksi[1][2]);
	new_axis.ksi[1][0] /= modul;
	new_axis.ksi[1][1] /= modul;
	new_axis.ksi[1][2] /= modul;

	// Third axis - vector producs
	new_axis.ksi[2][0] = new_axis.ksi[0][1] * new_axis.ksi[1][2] - new_axis.ksi[0][2] * new_axis.ksi[1][1];
	new_axis.ksi[2][1] = new_axis.ksi[0][2] * new_axis.ksi[1][0] - new_axis.ksi[0][0] * new_axis.ksi[1][2];
	new_axis.ksi[2][2] = new_axis.ksi[0][0] * new_axis.ksi[1][1] - new_axis.ksi[0][1] * new_axis.ksi[1][0];

	// Rotation matrix - random teta around normal
	float rot[3][3];

	rot[0][0] = cos(teta) + x * x * (1 - cos(teta));
	rot[0][1] = y * x * (1 - cos(teta)) + z * sin(teta);
	rot[0][2] = z * x * (1 - cos(teta)) - y * sin(teta);

	rot[1][0] = y * x * (1 - cos(teta)) - z * sin(teta);
	rot[1][1] = cos(teta) + y * y * (1 - cos(teta));
	rot[1][2] = y * z * (1 - cos(teta)) + x * sin(teta);

	rot[2][0] = x * z * (1 - cos(teta)) + y * sin(teta);
	rot[2][1] = y * z * (1 - cos(teta)) - x * sin(teta);
	rot[2][2] = cos(teta) + z * z * (1 - cos(teta));

	// Rotate
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			random_axis[node_num].ksi[i][j] = 0;
			for(int k = 0; k < 3; k++) {
				random_axis[node_num].ksi[i][j] += rot[k][j] * new_axis.ksi[i][k];
			}
		}
	}

};

void GCM_Tetr_Plastic_Interpolation_1stOrder_Rotate_Axis::find_inverse_matrix(int node_num)
{
	// Find inverse matrix, use the fact that it is rotation matrix

	random_axis_inv[node_num].ksi[0][0] = random_axis[node_num].ksi[0][0];
	random_axis_inv[node_num].ksi[1][0] = random_axis[node_num].ksi[0][1];
	random_axis_inv[node_num].ksi[2][0] = random_axis[node_num].ksi[0][2];

	random_axis_inv[node_num].ksi[0][1] = random_axis[node_num].ksi[1][0];
	random_axis_inv[node_num].ksi[1][1] = random_axis[node_num].ksi[1][1];
	random_axis_inv[node_num].ksi[2][1] = random_axis[node_num].ksi[1][2];

	random_axis_inv[node_num].ksi[0][2] = random_axis[node_num].ksi[2][0];
	random_axis_inv[node_num].ksi[1][2] = random_axis[node_num].ksi[2][1];
	random_axis_inv[node_num].ksi[2][2] = random_axis[node_num].ksi[2][2];

};

