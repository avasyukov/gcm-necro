#include <iomanip>

#include "TetrMesh_1stOrder.h"

TetrMesh_1stOrder::TetrMesh_1stOrder()
{
	mesh_type.assign("Tetrahedron mesh 1st order");
	T = gsl_matrix_alloc(3, 3);
	S = gsl_matrix_alloc(3, 3);
	P = gsl_permutation_alloc(3);
	mesh_min_h = -1;
};

TetrMesh_1stOrder::~TetrMesh_1stOrder()
{
	gsl_matrix_free(T);
	gsl_matrix_free(S);
	gsl_permutation_free(P);
	clear_data();
};

void TetrMesh_1stOrder::clear_data()
{
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].contact_data != NULL)
			free(nodes[i].contact_data);
		if(nodes[i].elements != NULL) {
			nodes[i].elements->clear();
			delete nodes[i].elements;
		}
		if(nodes[i].border_elements != NULL) {
			nodes[i].border_elements->clear();
			delete nodes[i].border_elements;
		}
	}
	nodes.clear();
	new_nodes.clear();
	tetrs.clear();
	border.clear();
	mesh_min_h = -1;
};

void TetrMesh_1stOrder::clear_contact_state()
{
	for(int i = 0; i < nodes.size(); i++)
		clear_contact_data(&nodes[i]);
};

//TODO - ugly hack - think about if we need it and how to implement it
void TetrMesh_1stOrder::add_node(ElasticNode* node)
{
	ElasticNode new_node = *node;
	new_node.mesh = this;
	nodes.push_back(new_node);
	new_nodes.push_back(new_node);
};

//TODO - ugly hack - think about if we need it and how to implement it
void TetrMesh_1stOrder::add_tetr(Tetrahedron_1st_order* tetr)
{
	Tetrahedron_1st_order new_tetr = *tetr;
	tetrs.push_back(new_tetr);
};

int TetrMesh_1stOrder::pre_process_mesh()
{
	*logger < "Preprocessing mesh started...";

//	*logger < "COORDS";
//	for(int i = 0; i < nodes.size(); i++)
//		*logger << i << " " << nodes[i].coords[0] << " " << nodes[i].coords[1] << " " < nodes[i].coords[2];
//	*logger < "IDX";
//	for(int i = 0; i < tetrs.size(); i++)
//		*logger << i << " " << tetrs[i].vert[0] << " " << tetrs[i].vert[1] << " " << tetrs[i].vert[2] << " " < tetrs[i].vert[3];
	
	// Just to ensure loaded border was dropped
	border.clear();

	// Guaranteed allowed step
	float step_h = get_max_h() / 4; // TODO avoid magick number


	*logger < "Checking numbering";

	// Check if internal numbers of nodes are the same as numbers in array
	// We need it in future to perform quick access to nodes in array
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].local_num != i)
			throw GCMException( GCMException::MESH_EXCEPTION, "Invalid node numbering");
	}
	for(int i = 0; i < tetrs.size(); i++)
	{
		if(tetrs[i].local_num != i)
			throw GCMException( GCMException::MESH_EXCEPTION, "Invalid tetrahedron numbering");
	}
	for(int i = 0; i < border.size(); i++)
	{
		if(border[i].local_num != i)
			throw GCMException( GCMException::MESH_EXCEPTION, "Invalid triangle numbering");
	}

	*logger < "Building volume reverse lookups";

	// Init vectors for "reverse lookups" of tetrahedrons current node is a member of.
	for(int i = 0; i < nodes.size(); i++) { nodes[i].elements = new vector<int>; }

	// Go through all the tetrahedrons
	for(int i = 0; i < tetrs.size(); i++)
	{
		// For all verticles
		for(int j = 0; j < 4; j++)
		{
			// Push to data of nodes the number of this tetrahedron
			nodes[tetrs[i].vert[j]].elements->push_back(i);
		}
	}

	*logger < "Looking for unused nodes";

	// Check all the nodes and find 'unused'

	// Case 1 - node has no connections at all
	for(int i = 0; i < nodes.size(); i++)
		if( (nodes[i].elements)->size() == 0 )
			nodes[i].setUsed (false);

	// Case 2 - remote ones that have connections only with remote ones
	for(int i = 0; i < nodes.size(); i++)
	{
		// If node is remote
		if(nodes[i].isRemote ())
		{
			int count = 0;
			// Check tetrahedrons it is a member of
			for(int j = 0; j < nodes[i].elements->size(); j++)
			{
				// Check verticles
				for(int k = 0; k < 4; k++)
				{
					// If it is local - count++
					if(nodes[tetrs[nodes[i].elements->at(j)].vert[k]].isLocal ()) { count++; }
				}
			}
			// If remote node is NOT connected with at least one local - it is unused one
			if(count == 0)
			{
				nodes[i].setUsed (false);
			}
		}
	}

	// Prepare border data

	float solid_angle;
	float ftmp;

	*logger < "Looking for border nodes using angles";

	// Check border using solid angle comparation with 4*PI
	for(int i = 0; i < nodes.size(); i++)
	{
		if( nodes[i].isLocal ())
		{
			solid_angle = 0;
			!nodes[i].isBorder ();//Why do we need it??
			for(int j = 0; j < nodes[i].elements->size(); j++)
			{
				ftmp = get_solid_angle(i, nodes[i].elements->at(j));
				if(ftmp < 0) { return -1; }
				solid_angle += ftmp;
			}
			if( (4 * qm_engine.PI - solid_angle) > qm_engine.PI / 100 ) { // TODO avoid magick number
				nodes[i].setIsBorder (true);
				nodes[i].contact_data = (contact_state*) malloc(sizeof(contact_state));
				clear_contact_data(&nodes[i]);
			}
		}
	}

	*logger < "Constructing border triangles";

	// Check all tetrs and construct border triangles
	for(int i = 0; i < tetrs.size(); i++)
	{
		check_triangle_to_be_border(tetrs[i].vert[0], tetrs[i].vert[1], tetrs[i].vert[2], tetrs[i].vert[3], step_h);
		check_triangle_to_be_border(tetrs[i].vert[0], tetrs[i].vert[1], tetrs[i].vert[3], tetrs[i].vert[2], step_h);
		check_triangle_to_be_border(tetrs[i].vert[0], tetrs[i].vert[2], tetrs[i].vert[3], tetrs[i].vert[1], step_h);
		check_triangle_to_be_border(tetrs[i].vert[1], tetrs[i].vert[2], tetrs[i].vert[3], tetrs[i].vert[0], step_h);
	}
	*logger < "Building surface reverse lookups";

	// Init vectors for "reverse lookups" of border triangles current node is a member of.
	for(int i = 0; i < nodes.size(); i++) { nodes[i].border_elements = new vector<int>; }

	// Go through all the triangles and push to data of nodes the number of this triangle
	for(int i = 0; i < border.size(); i++)
		for(int j = 0; j < 3; j++)
			nodes[border[i].vert[j]].border_elements->push_back(i);

	*logger < "Checking nodes outer normals";

	// Normal vector
	float normal[3];
	// Displacement
	float dx[3];

	// Check all nodes
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].isBorder ()) {

			find_border_node_normal(i, &normal[0], &normal[1], &normal[2], step_h);

			// Displacement along normal
			dx[0] = step_h * normal[0];
			dx[1] = step_h * normal[1];
			dx[2] = step_h * normal[2];

			// Check if we are inside of the body moving towards normal
			if( find_owner_tetr(&nodes[i], -dx[0], -dx[1], -dx[2]) == NULL ) {
				// Smth bad happens
//				*logger << "Can not find outer normal\n";
//				*logger << nodes[i].coords[0] << " " << nodes[i].coords[1] << " " < nodes[i].coords[2];
//				*logger << normal[0] << " " << normal[1] << " " < normal[2];
//				*logger << dx[0] << " " << dx[1] << " " < dx[2];
//				throw GCMException( GCMException::MESH_EXCEPTION, "Can not find outer normal");
			} else {
				// Outside body - normal is outer - do nothing
			}
		}

	}

	// TODO - scale, rotate, translate, etc - after it is made configurable via xml

	if (nodes.size())
	{
		*logger < "Creating outline";

		// Create outline
		for(int j = 0; j < 3; j++)
			outline.min_coords[j] = outline.max_coords[j] = nodes[0].coords[j];

		for(int i = 0; i < nodes.size(); i++) {
			for(int j = 0; j < 3; j++) {
				if(nodes[i].coords[j] > outline.max_coords[j])
					outline.max_coords[j] = nodes[i].coords[j];
				if(nodes[i].coords[j] < outline.min_coords[j])
					outline.min_coords[j] = nodes[i].coords[j];
			}
		}
	}

	*logger < "Preprocessing mesh done.";

	return 0;
};

int TetrMesh_1stOrder::check_triangle_to_be_border(int vert1, int vert2, int vert3, int tetr_vert, float step_h)
{
	int v1 = vert1;
	int v2 = vert2;
	int v3 = vert3;

	if( (!nodes[v1].isBorder ()) || (!nodes[v2].isBorder ()) || (!nodes[v3].isBorder ())
		|| (!nodes[v1].isLocal () ) || (!nodes[v2].isLocal ()) || (!nodes[v3].isLocal ()) )
		return -1;

	Triangle new_triangle;

	// Normal vector
	float normal[3];
	// Displacement
	float dx[3];

	float tri_center[3];
	float tetr_vert_direction[3];
	float plane_move[3];
	// float plane_move_mod;

	find_border_elem_normal(v1, v2, v3, &normal[0], &normal[1], &normal[2]);

	// Triangle center
	for(int i = 0; i < 3; i++)
		tri_center[i] = (nodes[v1].coords[i] + nodes[v2].coords[i] + nodes[v3].coords[i]) / 3;

	// Direction from triangle center to the verticle of tetrahedron
	for(int i = 0; i < 3; i++)
		tetr_vert_direction[i] = nodes[tetr_vert].coords[i] - tri_center[i];

	// Check if normal is co-linear with tetr verticle direction
	if(normal[0] * tetr_vert_direction[0] + normal[1] * tetr_vert_direction[1] + normal[2] * tetr_vert_direction[2] > 0)
	{
		// If they are - invert normal because only opposite direction can be outer
		for(int i = 0; i < 3; i++)
			normal[i] = -normal[i];
		// And swap verticles order to match new direction
		v2 = vert3;
		v3 = vert2;
	}

	// Displacement along potential outer normal
	for(int i = 0; i < 3; i++)
		dx[i] = step_h * normal[i];

	// Move from v1 to triangle center
	for(int i = 0; i < 3; i++)
		plane_move[i] = tri_center[i] - nodes[v1].coords[i];

	// Check if we are outside of the body moving from triangle center along normal ...
	if( find_owner_tetr(&nodes[v1], dx[0] + plane_move[0], dx[1] + plane_move[1], dx[2] + plane_move[2]) == NULL )
	{
		// ... add trianlge to border
		new_triangle.local_num = border.size();
		new_triangle.vert[0] = v1; new_triangle.vert[1] = v2; new_triangle.vert[2] = v3;
		border.push_back(new_triangle);
		return 0;
	} else {
		return -1;
	}

	return 0;
};

int TetrMesh_1stOrder::find_border_elem_normal(int border_element_index, float* x, float* y, float* z)
{
	int i = border_element_index;
	return find_border_elem_normal(border[i].vert[0], border[i].vert[1], border[i].vert[2], x, y, z);
};

int TetrMesh_1stOrder::find_border_elem_normal(float *p1, float *p2, float *p3, float *x, float *y, float *z)
{
	// Normal vector
	float normal[3];

	// Tmp vectors
	float v[2][3];

	// Tmp length
	float l;

	// Vector from vert '0' to vert '1'
	v[0][0] = p2[0] - p1[0];
	v[0][1] = p2[1] - p1[1];
	v[0][2] = p2[2] - p1[2];

	// Vector from vert '0' to vert '2'
	v[1][0] = p3[0] - p1[0];
	v[1][1] = p3[1] - p1[1];
	v[1][2] = p3[2] - p1[2];

	// Normal calculated as vector product
	normal[0] = v[0][1] * v[1][2] - v[0][2] * v[1][1];
	normal[1] = v[0][2] * v[1][0] - v[0][0] * v[1][2];
	normal[2] = v[0][0] * v[1][1] - v[0][1] * v[1][0];

	// Normal length is 1
	l = sqrt( normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2] );
	normal[0] /= l;
	normal[1] /= l;
	normal[2] /= l;

	*x = normal[0];
	*y = normal[1];
	*z = normal[2];

	return 0;
}

int TetrMesh_1stOrder::find_border_elem_normal(int v1, int v2, int v3, float* x, float* y, float* z)
{
	return find_border_elem_normal(nodes[v1].coords, nodes[v2].coords, nodes[v3].coords, x, y, z);
};

int TetrMesh_1stOrder::find_border_node_normal(int border_node_index, float* x, float* y, float* z)
{
	return find_border_node_normal(border_node_index, x, y, z, get_min_h());
};

int TetrMesh_1stOrder::find_border_node_normal(int border_node_index, float* x, float* y, float* z, float min_h /* Remove h from parameters */)
{
	float final_normal[3];
	final_normal[0] = 0;
	final_normal[1] = 0;
	final_normal[2] = 0;

	float cur_normal[3];

	int count = nodes[border_node_index].border_elements->size();

	for(int i = 0; i < count; i++) {
		find_border_elem_normal( (int) nodes[border_node_index].border_elements->at(i), &cur_normal[0], &cur_normal[1], &cur_normal[2] );
		final_normal[0] += cur_normal[0];	final_normal[1] += cur_normal[1];	final_normal[2] += cur_normal[2];
	}

	final_normal[0] /= count;
	final_normal[1] /= count;
	final_normal[2] /= count;

	float l = sqrt( final_normal[0] * final_normal[0] + final_normal[1] * final_normal[1] + final_normal[2] * final_normal[2] );
	final_normal[0] /= l;
        final_normal[1] /= l;
        final_normal[2] /= l;

	*x = final_normal[0];
	*y = final_normal[1];
	*z = final_normal[2];
	
       if( find_owner_tetr(&nodes[border_node_index], - min_h * final_normal[0], - min_h * final_normal[1], - min_h * final_normal[2]) != NULL )
               return 0;

       *logger < "DEBUG: using find_border_node_normal last chance option";
       int tetr_num = nodes[border_node_index].elements->at(0);

       float vector_end[3];
       for(int i = 0; i < 3; i++) {
               vector_end[i] = 0;
               for(int j = 0; j < 4; j++)
                       vector_end[i] += nodes[ tetrs[tetr_num].vert[j] ].coords[i];
               vector_end[i] /= 4;
               final_normal[i] = nodes[border_node_index].coords[i] - vector_end[i];
       }

       l = sqrt( final_normal[0] * final_normal[0] + final_normal[1] * final_normal[1] + final_normal[2] * final_normal[2] );
       final_normal[0] /= l;
        final_normal[1] /= l;
        final_normal[2] /= l;

       *x = final_normal[0];
       *y = final_normal[1];
       *z = final_normal[2];

	

	// If we are inside the body moving against normal - everything ok, just return
	if( find_owner_tetr(&nodes[border_node_index], - min_h * final_normal[0], - min_h * final_normal[1], - min_h * final_normal[2]) != NULL )
		return 0;

	// If we are here - it means the node is on some 'sharp' edge
	// In this case normal calculated as average of elements' normals is out of the body
	// Failback option - do not use average from all elements, get average between min-max only
	*logger < "DEBUG: using find_border_node_normal failback option";

	float min_coords[3];
	float max_coords[3];

	find_border_elem_normal( (int) nodes[border_node_index].border_elements->at(0), &cur_normal[0], &cur_normal[1], &cur_normal[2] );

	for(int i = 0; i < 3; i++) {
		min_coords[i] = cur_normal[i];		max_coords[i] = cur_normal[i];
	}

	for(int i = 1; i < count; i++) {
		find_border_elem_normal( (int) nodes[border_node_index].border_elements->at(i), &cur_normal[0], &cur_normal[1], &cur_normal[2] );
		for(int i = 0; i < 3; i++) {
			if( max_coords[i] < cur_normal[i] )
				max_coords[i] = cur_normal[i];
			if( min_coords[i] > cur_normal[i] )
				min_coords[i] = cur_normal[i];
		}
	}

	for(int i = 0; i < 3; i++)
		final_normal[i] = min_coords[i] + max_coords[i];

	l = sqrt( final_normal[0] * final_normal[0] + final_normal[1] * final_normal[1] + final_normal[2] * final_normal[2] );
	final_normal[0] /= l;
        final_normal[1] /= l;
        final_normal[2] /= l;

	*x = final_normal[0];
	*y = final_normal[1];
	*z = final_normal[2];

	if( find_owner_tetr(&nodes[border_node_index], - min_h * final_normal[0], - min_h * final_normal[1], - min_h * final_normal[2]) == NULL )
		throw GCMException( GCMException::MESH_EXCEPTION, "Can not create reasonable normal for node - is the border too sharp?");
}

// TODO move actual file and string operations into TaskPreparator or MshFileReader
int TetrMesh_1stOrder::load_ani3d_out_file(char* ani3d_out_file_name)
{
	int tmp_int;

	int number_of_nodes;
	int number_of_elements;
	ElasticNode new_node;
	Tetrahedron_1st_order new_tetr;

	ifstream ani3d_out_infile;

	ani3d_out_infile.open(ani3d_out_file_name, ifstream::in);
	if(!ani3d_out_infile.is_open())
		throw GCMException( GCMException::MESH_EXCEPTION, "Can not open ani3d_out file");

	*logger < "INFO: TetrMesh_1stOrder::load_ani3d_out_files - Reading file...";

	ani3d_out_infile >> number_of_nodes;

	for(int i = 0; i < number_of_nodes; i++)
	{
		// Zero all values
		new_node.local_num = new_node.remote_num = new_node.absolute_num = -1;
		new_node.local_zone_num = new_node.remote_zone_num = -1;
		new_node.coords[0] = new_node.coords[1] = new_node.coords[2] = 0;
		new_node.fixed_coords[0] = new_node.fixed_coords[1] = new_node.fixed_coords[2] = 0;
		new_node.la = new_node.mu = new_node.rho = 0;
		new_node.values[0] = new_node.values[1] = new_node.values[2] = 0;
		new_node.values[3] = new_node.values[4] = new_node.values[5] = 0;
		new_node.values[6] = new_node.values[7] = new_node.values[8] = 0;
		new_node.elements = NULL;
		new_node.contact_data = NULL;
		new_node.local_basis = NULL;
		new_node.setIsBorder (false);
		new_node.addOwner (GCM);
		new_node.setContactType (Free);
		new_node.setUsed (true);
		new_node.volume_calculator = NULL;
		new_node.border_condition = NULL;
		new_node.contact_condition = NULL;
		for(int z = 0; z < 8; z++)
			new_node.destruction_criterias[z] = 0;

		new_node.local_num = i;
		new_node.absolute_num = new_node.local_num;
		new_node.local_zone_num = zone_num;
		ani3d_out_infile >> new_node.coords[0] >> new_node.coords[1] >> new_node.coords[2];
		new_node.setPlacement (Local);

		new_node.mesh = this;
		nodes.push_back(new_node);
		new_nodes.push_back(new_node);
	}

	ani3d_out_infile >> number_of_elements;

	for(int i = 0; i < number_of_elements; i++)
	{
		ani3d_out_infile >> new_tetr.vert[0] >> new_tetr.vert[1] >> new_tetr.vert[2] >> new_tetr.vert[3] >>tmp_int;

		if( (new_tetr.vert[0] <= 0) || (new_tetr.vert[1] <= 0) || (new_tetr.vert[2] <= 0) || (new_tetr.vert[3] <= 0) )
			throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

		new_tetr.vert[0]--; new_tetr.vert[1]--; new_tetr.vert[2]--; new_tetr.vert[3]--;
		new_tetr.local_num = i;

		tetrs.push_back(new_tetr);
	}

	*logger < "INFO: TetrMesh_1stOrder::load_ani3d_out_files - File read.";

	ani3d_out_infile.close();

	return 0;
};

// TODO move actual file and string operations into TaskPreparator or MshFileReader
int TetrMesh_1stOrder::load_node_ele_files(char* node_file_name, char* ele_file_name)
{
	int tmp_int;

	int number_of_nodes;
	int number_of_elements;
	ElasticNode new_node;
	Tetrahedron_1st_order new_tetr;

	ifstream node_infile;
	ifstream ele_infile;

	node_infile.open(node_file_name, ifstream::in);
	if(!node_infile.is_open())
		throw GCMException( GCMException::MESH_EXCEPTION, "Can not open node file");

	ele_infile.open(ele_file_name, ifstream::in);
	if(!ele_infile.is_open())
		throw GCMException( GCMException::MESH_EXCEPTION, "Can not open ele file");

	*logger < "INFO: TetrMesh_1stOrder::load_node_ele_files - Reading file...";

	node_infile >> number_of_nodes >> tmp_int >> tmp_int >> tmp_int;

	for(int i = 0; i < number_of_nodes; i++)
	{
		// Zero all values
		new_node.local_num = new_node.remote_num = new_node.absolute_num = -1;
		new_node.local_zone_num = new_node.remote_zone_num = -1;
		new_node.coords[0] = new_node.coords[1] = new_node.coords[2] = 0;
		new_node.fixed_coords[0] = new_node.fixed_coords[1] = new_node.fixed_coords[2] = 0;
		new_node.la = new_node.mu = new_node.rho = 0;
		new_node.values[0] = new_node.values[1] = new_node.values[2] = 0;
		new_node.values[3] = new_node.values[4] = new_node.values[5] = 0;
		new_node.values[6] = new_node.values[7] = new_node.values[8] = 0;
		new_node.elements = NULL;
		new_node.contact_data = NULL;
		new_node.local_basis = NULL;
		new_node.setIsBorder (false);
		new_node.addOwner (GCM);
		new_node.setContactType (Free);

		node_infile >> new_node.local_num;
		if(new_node.local_num > 0)
		{
			new_node.local_num--;
			new_node.absolute_num = new_node.local_num;
			node_infile >> new_node.coords[0] >> new_node.coords[1] >> new_node.coords[2];
			new_node.setPlacement (Local);
			new_node.setIsBorder (false);
			new_node.addOwner (GCM);
			new_node.setContactType (Free);
			new_node.contact_data = NULL;
			new_node.local_basis = NULL;
			// TODO set other values
		}
		else
		{
			throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");
		}
		new_node.mesh = this;
		nodes.push_back(new_node);
		new_nodes.push_back(new_node);
	}

	ele_infile >> number_of_elements >> tmp_int >> tmp_int;

	for(int i = 0; i < number_of_elements; i++)
	{
		ele_infile >> new_tetr.local_num >> new_tetr.vert[0] >> new_tetr.vert[1] >> new_tetr.vert[2] >> new_tetr.vert[3];

		if( (new_tetr.vert[0] <= 0) || (new_tetr.vert[1] <= 0) || (new_tetr.vert[2] <= 0) || (new_tetr.vert[3] <= 0) )
			throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

		new_tetr.vert[0]--; new_tetr.vert[1]--; new_tetr.vert[2]--; new_tetr.vert[3]--;
		new_tetr.local_num--;

		tetrs.push_back(new_tetr);
	}

	*logger < "INFO: TetrMesh_1stOrder::load_node_ele_files - File read.";

	node_infile.close();
	ele_infile.close();

	return 0;
};

// TODO move actual file and string operations into TaskPreparator or MshFileReader
int TetrMesh_1stOrder::load_gmv_file(char* file_name)
{
	string str;
	int tmp_int;
	int count;

	int number_of_nodes;
	int number_of_elements;
	ElasticNode new_node;
	Tetrahedron_1st_order new_tetr;

	ifstream infile;

	infile.open(file_name, ifstream::in);
	if(!infile.is_open())
		throw GCMException( GCMException::MESH_EXCEPTION, "Can not open gmv file");

	*logger < "INFO: TetrMesh_1stOrder::load_gmv_file - Reading file...";

	getline(infile, str);
	if(strcmp(str.c_str(),"gmvinput ascii") != 0)
		throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

	getline(infile, str);

	infile >> str >> number_of_nodes;

	count = 1;
	for(int i = 0; i < number_of_nodes; i++)
	{
		// Zero all values
		new_node.local_num = new_node.remote_num = new_node.absolute_num = -1;
		new_node.local_zone_num = new_node.remote_zone_num = -1;
		new_node.coords[0] = new_node.coords[1] = new_node.coords[2] = 0;
		new_node.fixed_coords[0] = new_node.fixed_coords[1] = new_node.fixed_coords[2] = 0;
		new_node.la = new_node.mu = new_node.rho = 0;
		new_node.values[0] = new_node.values[1] = new_node.values[2] = 0;
		new_node.values[3] = new_node.values[4] = new_node.values[5] = 0;
		new_node.values[6] = new_node.values[7] = new_node.values[8] = 0;
		new_node.elements = NULL;
		new_node.contact_data = NULL;
		new_node.local_basis = NULL;
		new_node.setIsBorder (false);
		new_node.addOwner (GCM);
		new_node.setContactType (Free);

		new_node.local_num = count - 1;
		new_node.absolute_num = new_node.local_num;
		infile >> new_node.coords[0] >> new_node.coords[1] >> new_node.coords[2];

		new_node.mesh = this;
		nodes.push_back(new_node);
		new_nodes.push_back(new_node);

		count++;
	}

	getline(infile, str);

	infile >> str >> number_of_elements;

	count = 1;
	for(int i = 0; i < number_of_elements; i++)
	{
		new_tetr.local_num = count;

		infile >> str >> tmp_int >> new_tetr.vert[0] >> new_tetr.vert[1] >> new_tetr.vert[2] >> new_tetr.vert[3];

		if( (new_tetr.vert[0] <= 0) || (new_tetr.vert[1] <= 0) || (new_tetr.vert[2] <= 0) || (new_tetr.vert[3] <= 0) )
			throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

		new_tetr.vert[0]--; new_tetr.vert[1]--; new_tetr.vert[2]--; new_tetr.vert[3]--;
		new_tetr.local_num--;

		tetrs.push_back(new_tetr);

		count++;
	}

	*logger < "INFO: TetrMesh_1stOrder::load_gmv_file - File read.";

	infile.close();

	return 0;
};

// TODO move actual file and string operations into TaskPreparator or MshFileReader
int TetrMesh_1stOrder::load_msh_file(char* file_name)
{
	string str;
	float tmp_float;
	int tmp_int;
	int number_of_nodes;
	int number_of_elements;
	ElasticNode new_node;
	Tetrahedron_1st_order new_tetr;
	Triangle new_triangle;

	ifstream infile;
	infile.open(file_name, ifstream::in);
	if(!infile.is_open())
		throw GCMException( GCMException::MESH_EXCEPTION, "Can not open msh file");

	*logger < "INFO: TetrMesh_1stOrder::load_msh_file - Reading file...";

	infile >> str;
	if(strcmp(str.c_str(),"$MeshFormat") != 0)
		throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

	infile >> tmp_float >> tmp_float >> tmp_float;

	infile >> str;
	if(strcmp(str.c_str(),"$EndMeshFormat") != 0)
		throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

	*logger < "INFO: TetrMesh_1stOrder::load_msh_file - Reading file. Header read.";

	infile >> str;
	if(strcmp(str.c_str(),"$Nodes") != 0)
		throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

	infile >> number_of_nodes;

	for(int i = 0; i < number_of_nodes; i++)
	{
		// Zero all values
		new_node.local_num = new_node.remote_num = new_node.absolute_num = -1;
		new_node.local_zone_num = new_node.remote_zone_num = -1;
		new_node.coords[0] = new_node.coords[1] = new_node.coords[2] = 0;
		new_node.fixed_coords[0] = new_node.fixed_coords[1] = new_node.fixed_coords[2] = 0;
		new_node.la = new_node.mu = new_node.rho = 0;
		new_node.values[0] = new_node.values[1] = new_node.values[2] = 0;
		new_node.values[3] = new_node.values[4] = new_node.values[5] = 0;
		new_node.values[6] = new_node.values[7] = new_node.values[8] = 0;
		new_node.elements = NULL;
		new_node.contact_data = NULL;
		new_node.local_basis = NULL;
		new_node.setIsBorder (false);
		new_node.addOwner (GCM);
		new_node.setContactType (Free);
		new_node.setUsed (true);
		new_node.volume_calculator = NULL;
		new_node.border_condition = NULL;
		new_node.contact_condition = NULL;
		for(int z = 0; z < 8; z++)
			new_node.destruction_criterias[z] = 0;

		infile >> new_node.local_num;
		if(new_node.local_num > 0)
		{
			new_node.local_num--;
			new_node.absolute_num = new_node.local_num;
			new_node.local_zone_num = zone_num;
			infile >> new_node.coords[0] >> new_node.coords[1] >> new_node.coords[2];
			new_node.setPlacement (Local);
			// TODO set other values
		}
		else if(new_node.local_num < 0)
		{
			new_node.local_num = -new_node.local_num;
			new_node.local_num--;
			new_node.absolute_num = new_node.local_num;
			new_node.local_zone_num = zone_num;
			infile >> new_node.remote_zone_num >> new_node.remote_num;
			new_node.remote_num--;
			new_node.setPlacement (Remote);
			// TODO set other values
		}
		else
		{
			throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");
		}
		new_node.mesh = this;
		nodes.push_back(new_node);
		new_nodes.push_back(new_node);
	}

	infile >> str;
	if(strcmp(str.c_str(),"$EndNodes") != 0)
		throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

	*logger < "INFO: TetrMesh_1stOrder::load_msh_file - Reading file. Nodes read.";

	infile >> str;
	if(strcmp(str.c_str(),"$Elements") != 0)
		throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

	infile >> number_of_elements;
	for(int i = 0; i < number_of_elements; i++)
	{
		infile >> tmp_int;
		infile >> tmp_int;
		if(tmp_int != 4) {
			getline(infile, str);
			continue;
		} else if (tmp_int == 4) {
			new_tetr.local_num = tetrs.size();
			infile >> tmp_int >> tmp_int >> tmp_int  
				>> new_tetr.vert[0] >> new_tetr.vert[1] >> new_tetr.vert[2] >> new_tetr.vert[3];

			if( (new_tetr.vert[0] <= 0) || (new_tetr.vert[1] <= 0) || (new_tetr.vert[2] <= 0) || (new_tetr.vert[3] <= 0) )
				throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

			new_tetr.vert[0]--; new_tetr.vert[1]--; new_tetr.vert[2]--; new_tetr.vert[3]--;

			tetrs.push_back(new_tetr);
		}
	}

	infile >> str;
	if(strcmp(str.c_str(),"$EndElements") != 0)
		throw GCMException( GCMException::MESH_EXCEPTION, "Wrong file format");

	*logger < "INFO: TetrMesh_1stOrder::load_msh_file - File read.";

	infile.close();

	return 0;
};

void TetrMesh_1stOrder::load_vtu_file(string filename)
{
	vtkXMLUnstructuredGridReader *xgr = vtkXMLUnstructuredGridReader::New();
	vtkUnstructuredGrid *g = vtkUnstructuredGrid::New();

	xgr->SetFileName(const_cast<char*>(filename.c_str()));
	xgr->Update();

	g = xgr->GetOutput();
	*logger << "Number of points: " < g->GetNumberOfPoints();
	*logger << "Number of cells: " < g->GetNumberOfCells();

	double v[3];
	vtkDoubleArray *vel = (vtkDoubleArray*) g->GetPointData()->GetArray("velocity");
	vel->SetNumberOfComponents(3);
	vtkDoubleArray *sxx = (vtkDoubleArray*) g->GetPointData()->GetArray("sxx");
	vtkDoubleArray *sxy = (vtkDoubleArray*) g->GetPointData()->GetArray("sxy");
	vtkDoubleArray *sxz = (vtkDoubleArray*) g->GetPointData()->GetArray("sxz");
	vtkDoubleArray *syy = (vtkDoubleArray*) g->GetPointData()->GetArray("syy");
	vtkDoubleArray *syz = (vtkDoubleArray*) g->GetPointData()->GetArray("syz");
	vtkDoubleArray *szz = (vtkDoubleArray*) g->GetPointData()->GetArray("szz");
	vtkDoubleArray *la = (vtkDoubleArray*) g->GetPointData()->GetArray("lambda");
	vtkDoubleArray *mu = (vtkDoubleArray*) g->GetPointData()->GetArray("mu");
	vtkDoubleArray *rho = (vtkDoubleArray*) g->GetPointData()->GetArray("rho");
	vtkDoubleArray *yield_limit = (vtkDoubleArray*) g->GetPointData()->GetArray("yieldLimit");
	vtkIntArray *flags = (vtkIntArray*) g->GetPointData()->GetArray("flags");

	ElasticNode new_node;

	for( int i = 0; i < g->GetNumberOfPoints(); i++ )
	{
		double* dp = g->GetPoint(i);

		new_node.local_num = new_node.absolute_num = i;
		new_node.remote_num = -1;
		new_node.local_zone_num = zone_num;
		new_node.remote_zone_num = -1;
		new_node.coords[0] = new_node.fixed_coords[0] = dp[0];
		new_node.coords[1] = new_node.fixed_coords[1] = dp[1];
		new_node.coords[2] = new_node.fixed_coords[2] = dp[2];
		new_node.la = la->GetValue(i);
		new_node.mu = mu->GetValue(i);
		new_node.rho = rho->GetValue(i);
		new_node.yield_limit = yield_limit->GetValue(i);
		vel->GetTupleValue(i, v);
		new_node.values[0] = v[0];
		new_node.values[1] = v[1];
		new_node.values[2] = v[2];
		new_node.values[3] = sxx->GetValue(i);
		new_node.values[4] = sxy->GetValue(i);
		new_node.values[5] = sxz->GetValue(i);
		new_node.values[6] = syy->GetValue(i);
		new_node.values[7] = syz->GetValue(i);
		new_node.values[8] = szz->GetValue(i);
		new_node.setFlags( flags->GetValue(i) );
		new_node.elements = NULL;
		new_node.contact_data = NULL;
		new_node.local_basis = NULL;
		new_node.volume_calculator = NULL;
		new_node.border_condition = NULL;
		new_node.contact_condition = NULL;
		for(int z = 0; z < 8; z++)
			new_node.destruction_criterias[z] = 0;
		new_node.mesh = this;
		nodes.push_back(new_node);
		new_nodes.push_back(new_node);
	}

	Tetrahedron_1st_order new_tetr;
	for( int i = 0; i < g->GetNumberOfCells(); i++ )
	{
		new_tetr.local_num = i;
		vtkTetra *vt = (vtkTetra*) g->GetCell(i);
		new_tetr.vert[0] = vt->GetPointId(0);
		new_tetr.vert[1] = vt->GetPointId(1);
		new_tetr.vert[2] = vt->GetPointId(2);
		new_tetr.vert[3] = vt->GetPointId(3);
		tetrs.push_back(new_tetr);
	}

//	throw GCMException( GCMException::MESH_EXCEPTION, "VTU reader not implemented yet");
};

float TetrMesh_1stOrder::tetr_h(int i)
{
		float min_h;
		float h;
		float area[4];
		// Find volume
		float vol = qm_engine.tetr_volume(
			nodes[tetrs[i].vert[1]].coords[0] - nodes[tetrs[i].vert[0]].coords[0], 
			nodes[tetrs[i].vert[1]].coords[1] - nodes[tetrs[i].vert[0]].coords[1],
			nodes[tetrs[i].vert[1]].coords[2] - nodes[tetrs[i].vert[0]].coords[2],
			nodes[tetrs[i].vert[2]].coords[0] - nodes[tetrs[i].vert[0]].coords[0],
			nodes[tetrs[i].vert[2]].coords[1] - nodes[tetrs[i].vert[0]].coords[1],
			nodes[tetrs[i].vert[2]].coords[2] - nodes[tetrs[i].vert[0]].coords[2],
			nodes[tetrs[i].vert[3]].coords[0] - nodes[tetrs[i].vert[0]].coords[0],
			nodes[tetrs[i].vert[3]].coords[1] - nodes[tetrs[i].vert[0]].coords[1],
			nodes[tetrs[i].vert[3]].coords[2] - nodes[tetrs[i].vert[0]].coords[2]
		);

		// Find area of first face (verticles - 0,1,2)
		area[0] = qm_engine.tri_area(
			nodes[tetrs[i].vert[1]].coords[0] - nodes[tetrs[i].vert[0]].coords[0],
			nodes[tetrs[i].vert[1]].coords[1] - nodes[tetrs[i].vert[0]].coords[1],
			nodes[tetrs[i].vert[1]].coords[2] - nodes[tetrs[i].vert[0]].coords[2],
                	nodes[tetrs[i].vert[2]].coords[0] - nodes[tetrs[i].vert[0]].coords[0],
			nodes[tetrs[i].vert[2]].coords[1] - nodes[tetrs[i].vert[0]].coords[1],
			nodes[tetrs[i].vert[2]].coords[2] - nodes[tetrs[i].vert[0]].coords[2]
		);

		// Find area of second face (verticles - 0,1,3)
		area[1] = qm_engine.tri_area(
			nodes[tetrs[i].vert[1]].coords[0] - nodes[tetrs[i].vert[0]].coords[0],
			nodes[tetrs[i].vert[1]].coords[1] - nodes[tetrs[i].vert[0]].coords[1],
			nodes[tetrs[i].vert[1]].coords[2] - nodes[tetrs[i].vert[0]].coords[2],
			nodes[tetrs[i].vert[3]].coords[0] - nodes[tetrs[i].vert[0]].coords[0],
			nodes[tetrs[i].vert[3]].coords[1] - nodes[tetrs[i].vert[0]].coords[1],
			nodes[tetrs[i].vert[3]].coords[2] - nodes[tetrs[i].vert[0]].coords[2]
		);

		// Find area of third face (verticles - 0,2,3)
		area[2] = qm_engine.tri_area(
			nodes[tetrs[i].vert[2]].coords[0] - nodes[tetrs[i].vert[0]].coords[0],
			nodes[tetrs[i].vert[2]].coords[1] - nodes[tetrs[i].vert[0]].coords[1],
			nodes[tetrs[i].vert[2]].coords[2] - nodes[tetrs[i].vert[0]].coords[2],
	                nodes[tetrs[i].vert[3]].coords[0] - nodes[tetrs[i].vert[0]].coords[0],
			nodes[tetrs[i].vert[3]].coords[1] - nodes[tetrs[i].vert[0]].coords[1],
			nodes[tetrs[i].vert[3]].coords[2] - nodes[tetrs[i].vert[0]].coords[2]
		);

		// Find area of third face (verticles - 1,2,3)
		area[3] = qm_engine.tri_area(
			nodes[tetrs[i].vert[2]].coords[0] - nodes[tetrs[i].vert[1]].coords[0],
			nodes[tetrs[i].vert[2]].coords[1] - nodes[tetrs[i].vert[1]].coords[1],
			nodes[tetrs[i].vert[2]].coords[2] - nodes[tetrs[i].vert[1]].coords[2],
                	nodes[tetrs[i].vert[3]].coords[0] - nodes[tetrs[i].vert[1]].coords[0],
			nodes[tetrs[i].vert[3]].coords[1] - nodes[tetrs[i].vert[1]].coords[1],
			nodes[tetrs[i].vert[3]].coords[2] - nodes[tetrs[i].vert[1]].coords[2]
		);

		// Check if all nodes are already loaded from other CPUs and tetrahadron is correct
		if(vol == 0){ *logger < tetrs[i].vert[0];*logger < tetrs[i].vert[1];*logger < tetrs[i].vert[2];*logger < tetrs[i].vert[3];
			throw GCMException( GCMException::MESH_EXCEPTION, "Tetr volume is zero");}

		for(int j = 0; j < 4; j++)
			if(area[j] == 0)
				throw GCMException( GCMException::MESH_EXCEPTION, "Tri area is zero");

		min_h = fabs(3*vol/area[0]);

		// Go through all faces of given tetrahedron
		for(int j = 1; j < 4; j++)
		{
			// Find height to this face
			h = fabs(3*vol/area[j]);
			// And check if we should update minimum height
			if(h < min_h) { min_h = h; }
		}

		return min_h;
};

// Finds minimum h over mesh
float TetrMesh_1stOrder::get_min_h()
{
	if(mesh_min_h < 0)
		calc_min_h();

	return mesh_min_h;
};

void TetrMesh_1stOrder::calc_min_h()
{
	float min_h = -1;
	float h;
	// Go through tetrahedrons
	for(int i = 0; i < tetrs.size(); i++)
	{
		if ( (!nodes[tetrs[i].vert[0]].isUsed ())
			|| (!nodes[tetrs[i].vert[1]].isUsed ())
			|| (!nodes[tetrs[i].vert[2]].isUsed ())
			|| (!nodes[tetrs[i].vert[3]].isUsed ()) )
			continue;

		// Get current h
		h = tetr_h(i);
		if(min_h < 0)
			min_h = h;
		// If it is negative - smth bad happens
		if(h < 0)
			throw GCMException( GCMException::MESH_EXCEPTION, "Min h is negative");

		// Otherwise - just find minimum
		if(h < min_h) { min_h = h; }
	}
	mesh_min_h = min_h;
};

// Finds maximum h over mesh
float TetrMesh_1stOrder::get_max_h()
{
	float h;
	float max_h = -1;
	// Go through tetrahedrons
	for(int i = 0; i < tetrs.size(); i++)
	{
		if ( (!nodes[tetrs[i].vert[0]].isUsed ())
			|| (!nodes[tetrs[i].vert[1]].isUsed ())
			|| (!nodes[tetrs[i].vert[2]].isUsed ())
			|| (!nodes[tetrs[i].vert[3]].isUsed ()) )
			continue;

		// Get current h
		h = tetr_h(i);
		// If it is negative - smth bad happens
		if(h < 0)
			throw GCMException( GCMException::MESH_EXCEPTION, "Max h is negative");

		// Otherwise - just find minimum
		if(h > max_h) { max_h = h; }
	}

	return max_h;
};

int TetrMesh_1stOrder::log_mesh_stats()
{
	float h;

	float max_h;
	float avg_h = 0;

	float hyst[10];
	hyst[0] = hyst[1] = hyst[2] = hyst[3] = hyst[4] = hyst[5] = hyst[6] = hyst[7] = hyst[8] = hyst[9] = 0;

	max_h = get_max_h();

	int num;

	*logger << "Number of nodes: " < nodes.size();
	*logger << "Number of tetrs: " < tetrs.size();

	for(int i = 0; i < tetrs.size(); i++)
	{
		if ( (!nodes[tetrs[i].vert[0]].isUsed ())
			|| (!nodes[tetrs[i].vert[1]].isUsed ())
			|| (!nodes[tetrs[i].vert[2]].isUsed ())
			|| (!nodes[tetrs[i].vert[3]].isUsed ()) )
			continue;

		// Get current h
		h = tetr_h(i);
		// If it is negative - smth bad happens
		if(h < 0)
			throw GCMException( GCMException::MESH_EXCEPTION, "Tetr h is negative");

		avg_h += h/tetrs.size();

		h = h / max_h;
		num = (int)(h/0.1);
		hyst[num]++;
	}

	*logger < "Mesh outline:";
	*logger << "MinX: " < outline.min_coords[0];
	*logger << "MaxX: " < outline.max_coords[0];
	*logger << "MinY: " < outline.min_coords[1];
	*logger << "MaxY: " < outline.max_coords[1];
	*logger << "MinZ: " < outline.min_coords[2];
	*logger << "MaxZ: " < outline.max_coords[2];

	*logger < "Mesh quality:";
	*logger << "Max H = " < get_max_h();
	*logger << "Min H = " < get_min_h();
	*logger << "Avg H = " < avg_h;
	*logger < "Histogramm:";
	for(int i = 0; i < 10; i++)
		*logger < hyst[i];

	return 0;
};

bool TetrMesh_1stOrder::point_in_tetr(float x, float y, float z, Tetrahedron_1st_order* tetr)
{
	return point_in_tetr(x, y, z, (Tetrahedron*) tetr);
};

// TODO - rewrite with calc_determ_pure_tetr and calc_determ_with_shift
bool TetrMesh_1stOrder::point_in_tetr(float x, float y, float z, Tetrahedron* tetr)
{
	#ifdef DEBUG_MESH_GEOMETRY
		*logger < "DEBUG: TetrMesh_1stOrder::point_in_tetr";
		*logger << "Point: x: " << x << " y: " << y << " z: " < z;
		*logger << "Tetr: " < tetr->local_num;
		for(int j = 0; j < 4; j++) {
			ElasticNode* tmp_node = get_node( tetr->vert[j] );
			*logger << "\t\tVert: " << j << " num: " << tmp_node->local_num << "\t"
				<< " x: " << tmp_node->coords[0]
				<< " y: " << tmp_node->coords[1]
				<< " z: " < tmp_node->coords[2];
		}
	}
	#endif

	float d1,d2;
	d1 = qm_engine.determinant(
		nodes[tetr->vert[1]].coords[0] - nodes[tetr->vert[0]].coords[0],
		nodes[tetr->vert[1]].coords[1] - nodes[tetr->vert[0]].coords[1],
		nodes[tetr->vert[1]].coords[2] - nodes[tetr->vert[0]].coords[2],
		nodes[tetr->vert[2]].coords[0] - nodes[tetr->vert[0]].coords[0],
		nodes[tetr->vert[2]].coords[1] - nodes[tetr->vert[0]].coords[1],
		nodes[tetr->vert[2]].coords[2] - nodes[tetr->vert[0]].coords[2],
		nodes[tetr->vert[3]].coords[0] - nodes[tetr->vert[0]].coords[0],
		nodes[tetr->vert[3]].coords[1] - nodes[tetr->vert[0]].coords[1],
		nodes[tetr->vert[3]].coords[2] - nodes[tetr->vert[0]].coords[2]
	);
	d2 = qm_engine.determinant(
		nodes[tetr->vert[1]].coords[0] - x,
		nodes[tetr->vert[1]].coords[1] - y,
		nodes[tetr->vert[1]].coords[2] - z,
		nodes[tetr->vert[2]].coords[0] - x,
		nodes[tetr->vert[2]].coords[1] - y,
		nodes[tetr->vert[2]].coords[2] - z,
		nodes[tetr->vert[3]].coords[0] - x,
		nodes[tetr->vert[3]].coords[1] - y,
		nodes[tetr->vert[3]].coords[2] - z
	);
	#ifdef DEBUG_MESH_GEOMETRY
	*logger << "\t\tStage1: d1: " << d1 << " d2: " < d2;
	#endif
	if(d1*d2 < 0) { return false; }

	d1 = qm_engine.determinant(
		nodes[tetr->vert[0]].coords[0] - nodes[tetr->vert[1]].coords[0],
		nodes[tetr->vert[0]].coords[1] - nodes[tetr->vert[1]].coords[1],
		nodes[tetr->vert[0]].coords[2] - nodes[tetr->vert[1]].coords[2],
		nodes[tetr->vert[2]].coords[0] - nodes[tetr->vert[1]].coords[0],
		nodes[tetr->vert[2]].coords[1] - nodes[tetr->vert[1]].coords[1],
		nodes[tetr->vert[2]].coords[2] - nodes[tetr->vert[1]].coords[2],
		nodes[tetr->vert[3]].coords[0] - nodes[tetr->vert[1]].coords[0],
		nodes[tetr->vert[3]].coords[1] - nodes[tetr->vert[1]].coords[1],
		nodes[tetr->vert[3]].coords[2] - nodes[tetr->vert[1]].coords[2]
	);
	d2 = qm_engine.determinant(
		nodes[tetr->vert[0]].coords[0] - x,
		nodes[tetr->vert[0]].coords[1] - y,
		nodes[tetr->vert[0]].coords[2] - z,
		nodes[tetr->vert[2]].coords[0] - x,
		nodes[tetr->vert[2]].coords[1] - y,
		nodes[tetr->vert[2]].coords[2] - z,
		nodes[tetr->vert[3]].coords[0] - x,
		nodes[tetr->vert[3]].coords[1] - y,
		nodes[tetr->vert[3]].coords[2] - z
	);
	#ifdef DEBUG_MESH_GEOMETRY
	*logger << "\t\tStage2: d1: " << d1 << " d2: " < d2;
	#endif
	if(d1*d2 < 0) { return false; }

	d1 = qm_engine.determinant(
		nodes[tetr->vert[0]].coords[0] - nodes[tetr->vert[2]].coords[0],
		nodes[tetr->vert[0]].coords[1] - nodes[tetr->vert[2]].coords[1],
		nodes[tetr->vert[0]].coords[2] - nodes[tetr->vert[2]].coords[2],
		nodes[tetr->vert[1]].coords[0] - nodes[tetr->vert[2]].coords[0],
		nodes[tetr->vert[1]].coords[1] - nodes[tetr->vert[2]].coords[1],
		nodes[tetr->vert[1]].coords[2] - nodes[tetr->vert[2]].coords[2],
		nodes[tetr->vert[3]].coords[0] - nodes[tetr->vert[2]].coords[0],
		nodes[tetr->vert[3]].coords[1] - nodes[tetr->vert[2]].coords[1],
		nodes[tetr->vert[3]].coords[2] - nodes[tetr->vert[2]].coords[2]
	);
	d2 = qm_engine.determinant(
		nodes[tetr->vert[0]].coords[0] - x,
		nodes[tetr->vert[0]].coords[1] - y,
		nodes[tetr->vert[0]].coords[2] - z,
		nodes[tetr->vert[1]].coords[0] - x,
		nodes[tetr->vert[1]].coords[1] - y,
		nodes[tetr->vert[1]].coords[2] - z,
		nodes[tetr->vert[3]].coords[0] - x,
		nodes[tetr->vert[3]].coords[1] - y,
		nodes[tetr->vert[3]].coords[2] - z
	);
	#ifdef DEBUG_MESH_GEOMETRY
	*logger << "\t\tStage3: d1: " << d1 << " d2: " < d2;
	#endif
	if(d1*d2 < 0) { return false; }

	d1 = qm_engine.determinant(
		nodes[tetr->vert[0]].coords[0] - nodes[tetr->vert[3]].coords[0],
		nodes[tetr->vert[0]].coords[1] - nodes[tetr->vert[3]].coords[1],
		nodes[tetr->vert[0]].coords[2] - nodes[tetr->vert[3]].coords[2],
		nodes[tetr->vert[1]].coords[0] - nodes[tetr->vert[3]].coords[0],
		nodes[tetr->vert[1]].coords[1] - nodes[tetr->vert[3]].coords[1],
		nodes[tetr->vert[1]].coords[2] - nodes[tetr->vert[3]].coords[2],
		nodes[tetr->vert[2]].coords[0] - nodes[tetr->vert[3]].coords[0],
		nodes[tetr->vert[2]].coords[1] - nodes[tetr->vert[3]].coords[1],
		nodes[tetr->vert[2]].coords[2] - nodes[tetr->vert[3]].coords[2]
	);
	d2 = qm_engine.determinant(
		nodes[tetr->vert[0]].coords[0] - x,
		nodes[tetr->vert[0]].coords[1] - y,
		nodes[tetr->vert[0]].coords[2] - z,
		nodes[tetr->vert[1]].coords[0] - x,
		nodes[tetr->vert[1]].coords[1] - y,
		nodes[tetr->vert[1]].coords[2] - z,
		nodes[tetr->vert[2]].coords[0] - x,
		nodes[tetr->vert[2]].coords[1] - y,
		nodes[tetr->vert[2]].coords[2] - z
	);
	#ifdef DEBUG_MESH_GEOMETRY
		*logger << "\t\tStage4: d1: " << d1 << " d2: " < d2;
	#endif
	if(d1*d2 < 0) { return false; }

	return true;
};

bool TetrMesh_1stOrder::point_in_tetr(int base_node_index, float dx, float dy, float dz, Tetrahedron_1st_order* tetr)
{
	return point_in_tetr(base_node_index, dx, dy, dz, (Tetrahedron*) tetr);
};

bool TetrMesh_1stOrder::point_in_tetr(int base_node_index, float dx, float dy, float dz, Tetrahedron* tetr)
{
	return point_in_tetr(base_node_index, dx, dy, dz, tetr, false);
};

bool TetrMesh_1stOrder::point_in_tetr(int base_node_index, float _dx, float _dy, float _dz, Tetrahedron* tetr, bool debug)
{

	// TODO FIXME - we use this to workaround very small displacements that fail due to floating point rounding errors
	float min_h = tetr_h(tetr->local_num);
	float delta = sqrt(_dx*_dx + _dy*_dy + _dz*_dz);

	float dx, dy, dz;

	if( (delta > 0) && (delta < min_h * 0.99) ) {
		dx = _dx * min_h * 0.99 / delta;
		dy = _dy * min_h * 0.99 / delta;
		dz = _dz * min_h * 0.99 / delta;
	} else {
		dx = _dx;
		dy = _dy;
		dz = _dz;
	}

	if(debug) {
		*logger << "MinH = " << min_h << " Delta = " < delta;
		*logger < "DEBUG: TetrMesh_1stOrder::point_in_tetr";
		*logger <<"Point: num: " << base_node_index << " _dx: " << _dx << " _dy: " << _dy << " _dz: " < _dz;
		*logger <<"Point: num: " << base_node_index << " dx: " << dx << " dy: " << dy << " dz: " < dz;
		*logger << "Tetr: " < tetr->local_num;
		for(int j = 0; j < 4; j++) {
			ElasticNode* tmp_node = get_node( tetr->vert[j] );
			*logger << "\t\tVert: " << j << " num: " << tmp_node->local_num << "\t"
				<< " x: " << tmp_node->coords[0]
				<< " y: " << tmp_node->coords[1]
				<< " z: " < tmp_node->coords[2];
		}
	}


	float Vol = fabs( qm_engine.tetr_volume(
		(nodes[tetr->vert[1]].coords[0])-(nodes[tetr->vert[0]].coords[0]),
		(nodes[tetr->vert[1]].coords[1])-(nodes[tetr->vert[0]].coords[1]),
		(nodes[tetr->vert[1]].coords[2])-(nodes[tetr->vert[0]].coords[2]),
		(nodes[tetr->vert[2]].coords[0])-(nodes[tetr->vert[0]].coords[0]),
		(nodes[tetr->vert[2]].coords[1])-(nodes[tetr->vert[0]].coords[1]),
		(nodes[tetr->vert[2]].coords[2])-(nodes[tetr->vert[0]].coords[2]),
		(nodes[tetr->vert[3]].coords[0])-(nodes[tetr->vert[0]].coords[0]),
		(nodes[tetr->vert[3]].coords[1])-(nodes[tetr->vert[0]].coords[1]),
		(nodes[tetr->vert[3]].coords[2])-(nodes[tetr->vert[0]].coords[2])
	) );

	float x = nodes[base_node_index].coords[0] + dx;
        float y = nodes[base_node_index].coords[1] + dy;
        float z = nodes[base_node_index].coords[2] + dz;

	float vols[4];

	vols[0] = fabs( qm_engine.tetr_volume(
		(nodes[tetr->vert[1]].coords[0])-x,
		(nodes[tetr->vert[1]].coords[1])-y,
		(nodes[tetr->vert[1]].coords[2])-z,
		(nodes[tetr->vert[2]].coords[0])-x,
		(nodes[tetr->vert[2]].coords[1])-y,
		(nodes[tetr->vert[2]].coords[2])-z,
		(nodes[tetr->vert[3]].coords[0])-x,
		(nodes[tetr->vert[3]].coords[1])-y,
		(nodes[tetr->vert[3]].coords[2])-z
	) );

	vols[1] = fabs( qm_engine.tetr_volume(
		(nodes[tetr->vert[0]].coords[0])-x,
		(nodes[tetr->vert[0]].coords[1])-y,
		(nodes[tetr->vert[0]].coords[2])-z,
		(nodes[tetr->vert[2]].coords[0])-x,
		(nodes[tetr->vert[2]].coords[1])-y,
		(nodes[tetr->vert[2]].coords[2])-z,
		(nodes[tetr->vert[3]].coords[0])-x,
		(nodes[tetr->vert[3]].coords[1])-y,
		(nodes[tetr->vert[3]].coords[2])-z
	) );

	vols[2] = fabs( qm_engine.tetr_volume(
		(nodes[tetr->vert[1]].coords[0])-x,
		(nodes[tetr->vert[1]].coords[1])-y,
		(nodes[tetr->vert[1]].coords[2])-z,
		(nodes[tetr->vert[0]].coords[0])-x,
		(nodes[tetr->vert[0]].coords[1])-y,
		(nodes[tetr->vert[0]].coords[2])-z,
		(nodes[tetr->vert[3]].coords[0])-x,
		(nodes[tetr->vert[3]].coords[1])-y,
		(nodes[tetr->vert[3]].coords[2])-z
	) );

	vols[3] = fabs( qm_engine.tetr_volume(
		(nodes[tetr->vert[1]].coords[0])-x,
		(nodes[tetr->vert[1]].coords[1])-y,
		(nodes[tetr->vert[1]].coords[2])-z,
		(nodes[tetr->vert[2]].coords[0])-x,
		(nodes[tetr->vert[2]].coords[1])-y,
		(nodes[tetr->vert[2]].coords[2])-z,
		(nodes[tetr->vert[0]].coords[0])-x,
		(nodes[tetr->vert[0]].coords[1])-y,
		(nodes[tetr->vert[0]].coords[2])-z
	) );

	if( vols[0] + vols[1] + vols[2] + vols[3] < Vol * 1.001 ) {
		if(debug) { *logger < "IN"; }
		return true;
	} else {
		if(debug) { *logger < "OUT"; }
		return false;
	}


	float d1,d2;


	if( triangleOrientationOk(tetr->vert[1], tetr->vert[2], tetr->vert[3]) ) {
		d1 = calc_determ_pure_tetr(tetr->vert[1], tetr->vert[2], tetr->vert[3], tetr->vert[0]);
		d2 = calc_determ_with_shift(tetr->vert[1], tetr->vert[2], tetr->vert[3], base_node_index, dx, dy, dz);
	} else {
		d1 = calc_determ_pure_tetr(tetr->vert[1], tetr->vert[3], tetr->vert[2], tetr->vert[0]);
		d2 = calc_determ_with_shift(tetr->vert[1], tetr->vert[3], tetr->vert[2], base_node_index, dx, dy, dz);
	}

	if(debug) 
		*logger << "\t\tStage1: d1: " << d1 << " d2: " < d2;

	if(d1*d2 < 0) { return false; }

	if( triangleOrientationOk(tetr->vert[0], tetr->vert[2], tetr->vert[3]) ) {
		d1 = calc_determ_pure_tetr(tetr->vert[0], tetr->vert[2], tetr->vert[3], tetr->vert[1]);
		d2 = calc_determ_with_shift(tetr->vert[0], tetr->vert[2], tetr->vert[3], base_node_index, dx, dy, dz);
	} else {
		d1 = calc_determ_pure_tetr(tetr->vert[0], tetr->vert[3], tetr->vert[2], tetr->vert[1]);
		d2 = calc_determ_with_shift(tetr->vert[0], tetr->vert[3], tetr->vert[2], base_node_index, dx, dy, dz);
	}

	if(debug) 
		*logger << "\t\tStage2: d1: " << d1 << " d2: " < d2;

	if(d1*d2 < 0) { return false; }

	if( triangleOrientationOk(tetr->vert[0], tetr->vert[1], tetr->vert[3]) ) {
		d1 = calc_determ_pure_tetr(tetr->vert[0], tetr->vert[1], tetr->vert[3], tetr->vert[2]);
		d2 = calc_determ_with_shift(tetr->vert[0], tetr->vert[1], tetr->vert[3], base_node_index, dx, dy, dz);
	} else {
		d1 = calc_determ_pure_tetr(tetr->vert[0], tetr->vert[3], tetr->vert[1], tetr->vert[2]);
		d2 = calc_determ_with_shift(tetr->vert[0], tetr->vert[3], tetr->vert[1], base_node_index, dx, dy, dz);
	}

	if(debug)
		*logger << "\t\tStage3: d1: " << d1 << " d2: " < d2;

	if(d1*d2 < 0) { return false; }

	if( triangleOrientationOk(tetr->vert[0], tetr->vert[1], tetr->vert[2]) ) {
		d1 = calc_determ_pure_tetr(tetr->vert[0], tetr->vert[1], tetr->vert[2], tetr->vert[3]);
		d2 = calc_determ_with_shift(tetr->vert[0], tetr->vert[1], tetr->vert[2], base_node_index, dx, dy, dz);
	} else {
		d1 = calc_determ_pure_tetr(tetr->vert[0], tetr->vert[2], tetr->vert[1], tetr->vert[3]);
		d2 = calc_determ_with_shift(tetr->vert[0], tetr->vert[2], tetr->vert[1], base_node_index, dx, dy, dz);
	}

	if(debug)
		*logger << "\t\tStage4: d1: " << d1 << " d2: " < d2;

	if(d1*d2 < 0) { return false; }

	return true;
};

Tetrahedron_1st_order* TetrMesh_1stOrder::find_owner_tetr(ElasticNode* node, float dx, float dy, float dz)
{
	return find_owner_tetr(node, dx, dy, dz, false);
};

Tetrahedron_1st_order* TetrMesh_1stOrder::find_owner_tetr(ElasticNode* node, float dx, float dy, float dz, bool debug)
{
	// TODO - clean the code - unnecessary 'x <-> dx' and node 'pointer <-> index' intermix
	int base_node = node->local_num;

	// We use this implementation now because tau is still min_h/max_L
	// In that case if not found in adjacent tetrs - not found at all and out of body
	// Checking adjacent tetrahedrons
//	for(int i = 0; i < (node->elements)->size(); i++)
//	{
//		if( point_in_tetr(base_node, dx, dy, dz, &tetrs[(node->elements)->at(i)], debug) )
//		{
//			return &tetrs[(node->elements)->at(i)];
//		}
//	}
//
//	return NULL;

	// FIXME TODO
	// The code below is about 10 times (!) slower compared with simple implementation above
	// Checked for the case tau = min_h/max_L

	float x = nodes[base_node].coords[0] + dx;
	float y = nodes[base_node].coords[1] + dy;
	float z = nodes[base_node].coords[2] + dz;

	// A square of distance between point in question and local node
	// Will be used to check if it is worth to continue search or point in question is out of body
	float R2 = dx * dx + dy * dy + dz * dz;

	#ifdef DEBUG_MESH_GEOMETRY
	*logger < "DEBUG: TetrMesh_1stOrder::find_owner_tetr";
	*logger << "\t\tR2: " < R2;
	#endif

	//TODO May be std::set will be better? It guarantees unique elements.
	vector<int> checked; // Already checked tetrahedrons
	vector<int> checking; // Tetrs we are checking just now
	vector<int> tmp_to_check; // tmp list, used to construct new 'checking' list

	checked.clear();
	checking.clear();
	tmp_to_check.clear();

	// If current tetrs are inside sphere of radius R or not. If not - we should stop search and return NULL
	bool inside_R = true;

	float local_x, local_y, local_z;

	for(int i = 0; i < (node->elements)->size(); i++)
		checking.push_back((node->elements)->at(i));

	while(inside_R)
	{
		inside_R = false;

		// Check current tetrs
		for(int i = 0; i < checking.size(); i++)
		{
			// If found - return result
			if( point_in_tetr(base_node, dx, dy, dz, &tetrs[checking[i]], debug) ) {
				#ifdef DEBUG_MESH_GEOMETRY
				*logger << "\tFound in tetr: " < checking[i];
				#endif

				return &tetrs[checking[i]];
			}

			// Check if this tetr is still inside sphere of radius R
			// If we have already found at least one tetr in the sphere - skip check
			if(!inside_R)
			{
				// For all verticles of current tetr
				for(int j = 0; j < 4; j++)
				{
					// Skip base node. Otherwise we'll get false positive insideR for the 1st and 2nd layers
					if( nodes[ tetrs[checking[i]].vert[j] ].local_num == base_node )
						break;

					local_x = nodes[ tetrs[checking[i]].vert[j] ].coords[0];
					local_y = nodes[ tetrs[checking[i]].vert[j] ].coords[1];
					local_z = nodes[ tetrs[checking[i]].vert[j] ].coords[2];
					// If its distance smaller than R
					if( ( (node->coords[0] - local_x) * (node->coords[0] - local_x) 
						+ (node->coords[1] - local_y) * (node->coords[1] - local_y) 
						+ (node->coords[2] - local_z) * (node->coords[2] - local_z) ) < R2 ) {
							inside_R = true;
					}
					// TODO FIXME In theory we whould check if sphere and tetr intersect.
					// Current check - if at least one vert is in sphere.
					// It can be wrong on turbo bad tetrs. It fails if
					// sphere goes 'through' tetr to next layer but does not include its verticles.
				}
			}
		}

		// If current layer is not in sphere - there is no need to create next layer - just return NULL
		if(!inside_R)
			return NULL;

		// If not found in current tetrs - create new lists

		// First - add checked tetrs to the 'checked' list
		for(int i = 0; i < checking.size(); i++)
			checked.push_back(checking[i]);

		// Create very long list of all tetrs adjacent to any tetr in 'checking' list
		// (with duplicates and tetrs already checked)
		tmp_to_check.clear();
		for(int i = 0; i < checking.size(); i++)
			for(int j = 0; j < 4; j++)
				for(int k = 0; k < nodes[ tetrs[checking[i]].vert[j] ].elements->size(); k++)
					tmp_to_check.push_back(nodes[ tetrs[checking[i]].vert[j] ].elements->at(k));

		// Remove duplicates
		sort( tmp_to_check.begin(), tmp_to_check.end() );
		tmp_to_check.erase( unique( tmp_to_check.begin(), tmp_to_check.end() ), tmp_to_check.end() );

		// Copy to 'checking' and removing elements that are already checked
		checking.clear();
		bool is_checked;
		for(int i = 0; i < tmp_to_check.size(); i++)
		{
			is_checked = false;
			for(int j = 0; j < checked.size(); j++)
				if(tmp_to_check[i] == checked[j])
					is_checked = true;
			if(!is_checked)
				checking.push_back(tmp_to_check[i]);
		}
	}

	return NULL;

};

Tetrahedron_1st_order* TetrMesh_1stOrder::find_border_cross(ElasticNode* node, float dx, float dy, float dz, ElasticNode* cross)
{
	// Not implemented yet
	Tetrahedron_1st_order* tetr = find_border_cross(node, dx, dy, dz, cross->coords);

	int sph_count = 0;

	if( tetr != NULL )
		for(int j = 0; j < 4; j++)
			if( nodes[ tetr->vert[j] ].isBorder() )
				if( nodes[ tetr->vert[j] ].isOwnedBy(SPH) )
					sph_count++;

	// FIXME - it is possible that we have 4 border nodes, isn't it?
	if( sph_count == 3 )
		cross->addOwner( SPH );

	return tetr;
};

Tetrahedron_1st_order* TetrMesh_1stOrder::find_border_cross(ElasticNode* node, float dx, float dy, float dz, float* coords)
{
	// TODO - clean the code - unnecessary 'x <-> dx' and node 'pointer <-> index' intermix
	int base_node = node->local_num;

	float x = nodes[base_node].coords[0] + dx;
	float y = nodes[base_node].coords[1] + dy;
	float z = nodes[base_node].coords[2] + dz;

	// A square of distance between point in question and local node
	// Will be used to check if it is worth to continue search or point in question is out of body
	float R2 = dx * dx + dy * dy + dz * dz;

	float dist = sqrt(R2);
	float direction[3];
	direction[0] = dx / dist;
	direction[1] = dy / dist;
	direction[2] = dz / dist;

	#ifdef DEBUG_MESH_GEOMETRY
	*logger < "DEBUG: TetrMesh_1stOrder::find_owner_tetr";
	*logger << "\t\tR2: " < R2;
	#endif

	//TODO May be std::set will be better? It guarantees unique elements.
	vector<int> checked; // Already checked tetrahedrons
	vector<int> checking; // Tetrs we are checking just now
	vector<int> tmp_to_check; // tmp list, used to construct new 'checking' list

	checked.clear();
	checking.clear();
	tmp_to_check.clear();

	// If current tetrs are inside sphere of radius R or not. If not - we should stop search and return NULL
	bool inside_R = true;

	float local_x, local_y, local_z;

	for(int i = 0; i < (node->elements)->size(); i++)
		checking.push_back((node->elements)->at(i));

	while(inside_R)
	{
		inside_R = false;

		// Check current tetrs
		for(int i = 0; i < checking.size(); i++)
		{
			int border_verts_count = 0;
			int border_verts[3];
			for(int j = 0; j < 4; j++)
				if( nodes[ tetrs[checking[i]].vert[j] ].isBorder ()) {
					border_verts[ border_verts_count ] = tetrs[checking[i]].vert[j];
					border_verts_count++;
				}
			// FIXME - it is possible that we have 4 border verts, isn't it?
			if( border_verts_count == 3 ) {

				*logger < "Border face found. Looking for intersection.";

				if( vector_intersects_triangle( 
						nodes[ border_verts[0] ].coords,
						nodes[ border_verts[1] ].coords,
						nodes[ border_verts[2] ].coords,
						node->coords,
						direction, dist, coords ) ) {
					*logger < "Intersection found.";
					return &tetrs[checking[i]];
				}
			}

			// Check if this tetr is still inside sphere of radius R
			// If we have already found at least one tetr in the sphere - skip check
			if(!inside_R)
			{
				// For all verticles of current tetr
				for(int j = 0; j < 4; j++)
				{
					// Skip base node. Otherwise we'll get false positive insideR for the 1st and 2nd layers
					if( nodes[ tetrs[checking[i]].vert[j] ].local_num == base_node )
						break;

					local_x = nodes[ tetrs[checking[i]].vert[j] ].coords[0];
					local_y = nodes[ tetrs[checking[i]].vert[j] ].coords[1];
					local_z = nodes[ tetrs[checking[i]].vert[j] ].coords[2];
					// If its distance smaller than R
					if( ( (node->coords[0] - local_x) * (node->coords[0] - local_x) 
						+ (node->coords[1] - local_y) * (node->coords[1] - local_y) 
						+ (node->coords[2] - local_z) * (node->coords[2] - local_z) ) < R2 ) {
							inside_R = true;
					}
					// TODO FIXME In theory we whould check if sphere and tetr intersect.
					// Current check - if at least one vert is in sphere.
					// It can be wrong on turbo bad tetrs. It fails if
					// sphere goes 'through' tetr to next layer but does not include its verticles.
				}
			}
		}

		// If current layer is not in sphere - there is no need to create next layer - just return NULL
		if(!inside_R)
			return NULL;

		// If not found in current tetrs - create new lists

		// First - add checked tetrs to the 'checked' list
		for(int i = 0; i < checking.size(); i++)
			checked.push_back(checking[i]);

		// Create very long list of all tetrs adjacent to any tetr in 'checking' list
		// (with duplicates and tetrs already checked)
		tmp_to_check.clear();
		for(int i = 0; i < checking.size(); i++)
			for(int j = 0; j < 4; j++)
				for(int k = 0; k < nodes[ tetrs[checking[i]].vert[j] ].elements->size(); k++)
					tmp_to_check.push_back(nodes[ tetrs[checking[i]].vert[j] ].elements->at(k));

		// Remove duplicates
		sort( tmp_to_check.begin(), tmp_to_check.end() );
		tmp_to_check.erase( unique( tmp_to_check.begin(), tmp_to_check.end() ), tmp_to_check.end() );

		// Copy to 'checking' and removing elements that are already checked
		checking.clear();
		bool is_checked;
		for(int i = 0; i < tmp_to_check.size(); i++)
		{
			is_checked = false;
			for(int j = 0; j < checked.size(); j++)
				if(tmp_to_check[i] == checked[j])
					is_checked = true;
			if(!is_checked)
				checking.push_back(tmp_to_check[i]);
		}
	}

	return NULL;

};


int TetrMesh_1stOrder::interpolate(ElasticNode* node, Tetrahedron* tetr)
{
	float Vol = qm_engine.tetr_volume(
		(nodes[tetr->vert[1]].coords[0])-(nodes[tetr->vert[0]].coords[0]),
		(nodes[tetr->vert[1]].coords[1])-(nodes[tetr->vert[0]].coords[1]),
		(nodes[tetr->vert[1]].coords[2])-(nodes[tetr->vert[0]].coords[2]),
		(nodes[tetr->vert[2]].coords[0])-(nodes[tetr->vert[0]].coords[0]),
		(nodes[tetr->vert[2]].coords[1])-(nodes[tetr->vert[0]].coords[1]),
		(nodes[tetr->vert[2]].coords[2])-(nodes[tetr->vert[0]].coords[2]),
		(nodes[tetr->vert[3]].coords[0])-(nodes[tetr->vert[0]].coords[0]),
		(nodes[tetr->vert[3]].coords[1])-(nodes[tetr->vert[0]].coords[1]),
		(nodes[tetr->vert[3]].coords[2])-(nodes[tetr->vert[0]].coords[2])
	);

	float factor[4];

	factor[0] = fabs( qm_engine.tetr_volume(
		(nodes[tetr->vert[1]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[1]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[1]].coords[2])-(node->coords[2]),
		(nodes[tetr->vert[2]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[2]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[2]].coords[2])-(node->coords[2]),
		(nodes[tetr->vert[3]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[3]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[3]].coords[2])-(node->coords[2])
	) / Vol);

	factor[1] = fabs( qm_engine.tetr_volume(
		(nodes[tetr->vert[0]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[0]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[0]].coords[2])-(node->coords[2]),
		(nodes[tetr->vert[2]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[2]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[2]].coords[2])-(node->coords[2]),
		(nodes[tetr->vert[3]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[3]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[3]].coords[2])-(node->coords[2])
	) / Vol);

	factor[2] = fabs( qm_engine.tetr_volume(
		(nodes[tetr->vert[1]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[1]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[1]].coords[2])-(node->coords[2]),
		(nodes[tetr->vert[0]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[0]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[0]].coords[2])-(node->coords[2]),
		(nodes[tetr->vert[3]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[3]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[3]].coords[2])-(node->coords[2])
	) / Vol);

	factor[3] = fabs( qm_engine.tetr_volume(
		(nodes[tetr->vert[1]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[1]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[1]].coords[2])-(node->coords[2]),
		(nodes[tetr->vert[2]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[2]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[2]].coords[2])-(node->coords[2]),
		(nodes[tetr->vert[0]].coords[0])-(node->coords[0]),
		(nodes[tetr->vert[0]].coords[1])-(node->coords[1]),
		(nodes[tetr->vert[0]].coords[2])-(node->coords[2])
	) / Vol);

	// If we see potential instability
	if(factor[0] + factor[1] + factor[2] + factor[3] > 1.0)
	{
		// If it is small - treat instability as minor and just 'smooth' it
		// TODO - think about it more carefully
		//if( point_in_tetr(node->local_num, node->coords[0], node->coords[1], node->coords[2], tetr) )
		if(factor[0] + factor[1] + factor[2] + factor[3] < 1.25)
		{
			float sum = factor[0] + factor[1] + factor[2] + factor[3];
			for(int i = 0; i < 4; i++)
				factor[i] = factor[i] * 0.995 / sum;
		}
		// If point is not in tetr - throw exception
		else
		{
			*logger << "\tfactor[0]=" << factor[0] << " factor[1]=" << factor[1] << " factor[2]=" << factor[2] 	<< " factor[3]=" << factor[3] << " Sum: " < factor[0] + factor[1] + factor[2] + factor[3];

			*logger << "\tnode.x[0]=" << node->coords[0] << " node.x[1]=" << node->coords[1] 
				<< " node.x[2]=" < node->coords[2];

			*logger << "\tv0.x[0]=" << nodes[tetr->vert[0]].coords[0] << " v0.x[1]=" << nodes[tetr->vert[0]].coords[1] << " v0.x[2]=" < nodes[tetr->vert[0]].coords[2];
						
			*logger << "\tv1.x[0]=" << nodes[tetr->vert[1]].coords[0] << " v1.x[1]=" << nodes[tetr->vert[1]].coords[1] << " v1.x[2]=" < nodes[tetr->vert[1]].coords[2];

			*logger << "\tv2.x[0]=" << nodes[tetr->vert[2]].coords[0] << " v2.x[1]=" << nodes[tetr->vert[2]].coords[1] << " v2.x[2]=" < nodes[tetr->vert[2]].coords[2];

			*logger << "\tv3.x[0]=" << nodes[tetr->vert[3]].coords[0] << " v3.x[1]=" << nodes[tetr->vert[3]].coords[1] << " v3.x[2]=" < nodes[tetr->vert[3]].coords[2];
			throw GCMException( GCMException::MESH_EXCEPTION, "Sum of factors is greater than 1.0");
		}
	}

	for (int i = 0; i < 3; i++)
	{
		node->fixed_coords[i] = node->coords[i];
	}

	for (int i = 0; i < 9; i++)
	{
		node->values[i] = (nodes[tetr->vert[0]].values[i]*factor[0] 
				+ nodes[tetr->vert[1]].values[i]*factor[1] 
				+ nodes[tetr->vert[2]].values[i]*factor[2] 
				+ nodes[tetr->vert[3]].values[i]*factor[3]);
	}

        node->la = (nodes[tetr->vert[0]].la*factor[0] + nodes[tetr->vert[1]].la*factor[1] 
			+ nodes[tetr->vert[2]].la*factor[2] + nodes[tetr->vert[3]].la*factor[3]);
        node->mu = (nodes[tetr->vert[0]].mu*factor[0] + nodes[tetr->vert[1]].mu*factor[1] 
			+ nodes[tetr->vert[2]].mu*factor[2] + nodes[tetr->vert[3]].mu*factor[3]);
        node->rho = (nodes[tetr->vert[0]].rho*factor[0] + nodes[tetr->vert[1]].rho*factor[1] 
			+ nodes[tetr->vert[2]].rho*factor[2] + nodes[tetr->vert[3]].rho*factor[3]);

	return 0;
};

float TetrMesh_1stOrder::get_max_possible_tau()
{
	float min_h = get_min_h();
	float max_l = 0;
	float l;

	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].isLocal ())
		{
			l = method->get_max_lambda(&nodes[i], this);
			if(l < 0) {
				throw GCMException( GCMException::MESH_EXCEPTION, "Got negative lambda");
			}
			if(l > max_l) { max_l = l; }
		}
	}
	return min_h/max_l;
};

int TetrMesh_1stOrder::do_next_part_step(float tau, int stage)
{
	// Border nodes
	*logger < "Processing border nodes";
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].isLocal ())
		{
			if( nodes[i].isBorder () )
				method->do_next_part_step(&nodes[i], &new_nodes[i], tau, stage, this);
		}
	}
	// Inner nodes
	*logger < "Processing inner nodes";
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].isLocal ())
		{
			if( nodes[i].isInner ())
				method->do_next_part_step(&nodes[i], &new_nodes[i], tau, stage, this);
		}
	}
	// Copy values
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].isLocal ())
		{
			for(int j = 0; j < 9; j++)
				nodes[i].values[j] = new_nodes[i].values[j];
		}
	}
	return 0;
};

void TetrMesh_1stOrder::move_coords(float tau)
{
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].isLocal ())
		{
			for(int j = 0; j < 3; j++)
			{
				// Move node
				nodes[i].coords[j] += nodes[i].values[j]*tau;
				// Move mesh outline if necessary
				if(nodes[i].coords[j] > outline.max_coords[j])
					outline.max_coords[j] = nodes[i].coords[j];
				if(nodes[i].coords[j] < outline.min_coords[j])
					outline.min_coords[j] = nodes[i].coords[j];
			}
		}
	}
	calc_min_h();
};

void TetrMesh_1stOrder::calc_destruction_criterias()
{
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].isLocal ())
		{
			nodes[i].max_compression = nodes[i].max_tension = nodes[i].max_shear = nodes[i].max_deviator = 0;

			// FIXME - we need main tensor components, not diagonal components
			if( nodes[i].values[3] > nodes[i].max_tension ) { nodes[i].max_tension = nodes[i].values[3]; }
			if( nodes[i].values[3] < nodes[i].max_compression ) { nodes[i].max_compression = nodes[i].values[3]; }
			if( nodes[i].values[6] > nodes[i].max_tension ) { nodes[i].max_tension = nodes[i].values[6]; }
			if( nodes[i].values[6] < nodes[i].max_compression ) { nodes[i].max_compression = nodes[i].values[6]; }
			if( nodes[i].values[8] > nodes[i].max_tension ) { nodes[i].max_tension = nodes[i].values[8]; }
			if( nodes[i].values[8] < nodes[i].max_compression ) { nodes[i].max_compression = nodes[i].values[8]; }
			nodes[i].max_compression = fabs( nodes[i].max_compression);

			if( fabs(nodes[i].values[4]) > nodes[i].max_shear ) { nodes[i].max_shear = fabs(nodes[i].values[4]); }
			if( fabs(nodes[i].values[5]) > nodes[i].max_shear ) { nodes[i].max_shear = fabs(nodes[i].values[5]); }
			if( fabs(nodes[i].values[7]) > nodes[i].max_shear ) { nodes[i].max_shear = fabs(nodes[i].values[7]); }

			nodes[i].max_deviator 
				= sqrt( ( (nodes[i].values[3] - nodes[i].values[6]) * (nodes[i].values[3] - nodes[i].values[6]) 
					+ (nodes[i].values[6] - nodes[i].values[8]) * (nodes[i].values[6] - nodes[i].values[8])
					+ (nodes[i].values[3] - nodes[i].values[8]) * (nodes[i].values[3] - nodes[i].values[8])
					+ 6 * ( (nodes[i].values[4]) * (nodes[i].values[4]) + (nodes[i].values[5]) * (nodes[i].values[5])
					+ (nodes[i].values[7]) * (nodes[i].values[7])) ) / 6 );

			if( nodes[i].max_tension > nodes[i].max_tension_history )
				nodes[i].max_tension_history = nodes[i].max_tension;
			if( nodes[i].max_compression > nodes[i].max_compression_history )
				nodes[i].max_compression_history = nodes[i].max_compression;
			if( nodes[i].max_shear > nodes[i].max_shear_history )
				nodes[i].max_shear_history = nodes[i].max_shear;
			if( nodes[i].max_deviator > nodes[i].max_deviator_history )
				nodes[i].max_deviator_history = nodes[i].max_deviator;
		}
	}
};

int TetrMesh_1stOrder::proceed_rheology()
{
	if(rheology == NULL)
		throw GCMException( GCMException::MESH_EXCEPTION, "No RC attached");

	// TODO if rheology is void we can skip this step
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].isLocal ())
		{
			rheology->do_calc(&nodes[i], &nodes[i]);
		}
	}
	return 0;
};

void TetrMesh_1stOrder::check_stresses(float tau)
{
	*logger << "Nodes: " < nodes.size();
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].isLocal ())
		{
// FIXME
//stresser->set_current_stress(&nodes[i], &nodes[i], tau);
			if( nodes[i].volume_calculator == NULL )
				nodes[i].volume_calculator = mesh_set->getDefaultVolumeCalculator();
			if( (nodes[i].border_condition == NULL) || ( ! nodes[i].border_condition->form->isActive(tau) ) )
				nodes[i].border_condition = mesh_set->getDefaultBorderCondition();
			if( (nodes[i].contact_condition == NULL) || ( ! nodes[i].contact_condition->form->isActive(tau) ) )
				nodes[i].contact_condition = mesh_set->getDefaultContactCondition();
		}
	}
};

int TetrMesh_1stOrder::run_mesh_filter()
{
	// Check all nodes
	for(int i = 0; i < nodes.size(); i++)
	{
		// Clear only local AND border (inner nodes have not caused problems till now)
		if( (nodes[i].isLocal ()) && (nodes[i].isBorder ()) )
		{
			// For each variable
			for(int j = 0; j < 9; j++)
			{
				bool alarm = true;
				float val = 0;
				int count = 0;
				// Check all tetrs around node
				for(int k = 0; k < nodes[i].elements->size(); k++) {
					// Check all verts
					for(int l = 0; l < 4; l++) {
						if( tetrs[nodes[i].elements->at(k)].vert[l] != nodes[i].local_num ) {
							// If value at this vert is the same sign - node is ok, stop check
							if( nodes[tetrs[nodes[i].elements->at(k)].vert[l]].values[j]
								 * nodes[i].values[j] >= 0 ) {
								 alarm = false;
							// If it is opposite sign - add it to possible 'smoothing' value
							} else {
								val += nodes[tetrs[nodes[i].elements->at(k)].vert[l]].values[j];
								count++;
							}
						}
					}
					if( !alarm )
						break;
				}
				// If all nodes around have different sign of this variable - just 'smooth' it in this node
				if( alarm ) {
					*logger < "INFO: TetrMesh_1stOrder::run_mesh_filter - node cleared";
					*logger << "\tNode " << nodes[i].local_num << ":"
							<< " x: " << nodes[i].coords[0]
							<< " y: " << nodes[i].coords[1]
							<< " z: " < nodes[i].coords[2];
					*logger << "\tVar: " << j 
							<< " Old value: " << nodes[i].values[j] 
							<< " New value: " < val / count;
					nodes[i].values[j] = val / count;
				}
			}
		}
	}
	return 0;
};

Tetrahedron* TetrMesh_1stOrder::get_tetrahedron(int index)
{
	return &tetrs[index];
};

ElasticNode* TetrMesh_1stOrder::get_node(int index)
{
	return &nodes[index];
};

float TetrMesh_1stOrder::calc_determ_pure_tetr(int node1, int node2, int node3, int ref_node)
{
	return qm_engine.determinant(
				nodes[node1].coords[0] - nodes[ref_node].coords[0],
				nodes[node1].coords[1] - nodes[ref_node].coords[1],
				nodes[node1].coords[2] - nodes[ref_node].coords[2],
				nodes[node2].coords[0] - nodes[ref_node].coords[0],
				nodes[node2].coords[1] - nodes[ref_node].coords[1],
				nodes[node2].coords[2] - nodes[ref_node].coords[2],
				nodes[node3].coords[0] - nodes[ref_node].coords[0],
				nodes[node3].coords[1] - nodes[ref_node].coords[1],
				nodes[node3].coords[2] - nodes[ref_node].coords[2] );
};

float TetrMesh_1stOrder::calc_determ_with_shift(int node1, int node2, int node3, int base_node, float dx, float dy, float dz)
{
	float x = nodes[base_node].coords[0] + dx;
	float y = nodes[base_node].coords[1] + dy;
	float z = nodes[base_node].coords[2] + dz;

	if(node1 == base_node) {
		return qm_engine.determinant(
				- dx,				- dy,				- dz,
				nodes[node2].coords[0] - x,	nodes[node2].coords[1] - y,	nodes[node2].coords[2] - z,
				nodes[node3].coords[0] - x,	nodes[node3].coords[1] - y,	nodes[node3].coords[2] - z );

	} else if (node2 == base_node) {
		return qm_engine.determinant(
				nodes[node1].coords[0] - x,	nodes[node1].coords[1] - y,	nodes[node1].coords[2] - z,
				- dx,				- dy,				- dz,
				nodes[node3].coords[0] - x,	nodes[node3].coords[1] - y,	nodes[node3].coords[2] - z );

	} else if (node3 == base_node) {
		return qm_engine.determinant(
				nodes[node1].coords[0] - x,	nodes[node1].coords[1] - y,	nodes[node1].coords[2] - z,
				nodes[node2].coords[0] - x,	nodes[node2].coords[1] - y,	nodes[node2].coords[2] - z,
				- dx,				- dy,				- dz			   );

	} else {
		return qm_engine.determinant(
				nodes[node1].coords[0] - x,	nodes[node1].coords[1] - y,	nodes[node1].coords[2] - z,
				nodes[node2].coords[0] - x,	nodes[node2].coords[1] - y,	nodes[node2].coords[2] - z,
				nodes[node3].coords[0] - x,	nodes[node3].coords[1] - y,	nodes[node3].coords[2] - z );
	}

};

bool TetrMesh_1stOrder::triangleOrientationOk(int node1, int node2, int node3)
{
	float x1 = nodes[node2].coords[0] - nodes[node1].coords[0];
	float x2 = nodes[node3].coords[0] - nodes[node1].coords[0];
	float y1 = nodes[node2].coords[1] - nodes[node1].coords[1];
	float y2 = nodes[node3].coords[1] - nodes[node1].coords[1];
	float z1 = nodes[node2].coords[2] - nodes[node1].coords[2];
	float z2 = nodes[node3].coords[2] - nodes[node1].coords[2];

	float xn = y1 * z2 - y2 * z1;
	if( xn > 0 ) { return true; }
	if( xn < 0 ) { return false; }

	float yn = -x1 * z2 + x2 * z1;
	if( yn > 0 ) { return true; }
	if( yn < 0 ) { return false; }

	float zn = x1 * y2 - x2 * y1;
	if( zn > 0 ) { return true; }
	if( zn < 0 ) { return false; }

	return false;
};

float TetrMesh_1stOrder::get_solid_angle(int node_index, int tetr_index)
{
	if( (tetrs[tetr_index].vert[0] != node_index) && (tetrs[tetr_index].vert[1] != node_index)
				&& (tetrs[tetr_index].vert[2] != node_index) && (tetrs[tetr_index].vert[3] != node_index) ) { return -1.0; }

	int vert[3];
	int count = 0;

	for(int i = 0; i < 4; i++)
		if(tetrs[tetr_index].vert[i] != node_index) {
			vert[count] = tetrs[tetr_index].vert[i];
			count++;
		}

	if(count != 3) { return -1.0; }

	return qm_engine.solid_angle(
			nodes[vert[0]].coords[0] - nodes[node_index].coords[0], 
			nodes[vert[0]].coords[1] - nodes[node_index].coords[1],
			nodes[vert[0]].coords[2] - nodes[node_index].coords[2],
			nodes[vert[1]].coords[0] - nodes[node_index].coords[0], 
			nodes[vert[1]].coords[1] - nodes[node_index].coords[1],
			nodes[vert[1]].coords[2] - nodes[node_index].coords[2],
			nodes[vert[2]].coords[0] - nodes[node_index].coords[0], 
			nodes[vert[2]].coords[1] - nodes[node_index].coords[1],
			nodes[vert[2]].coords[2] - nodes[node_index].coords[2]
	);
};

// For algo - see http://www.cs.princeton.edu/courses/archive/fall00/cs426/lectures/raycast/sld016.htm and use the brain
bool TetrMesh_1stOrder::vector_intersects_triangle(float *p1, float *p2, float *p3, float *p0, float *v, float l, float *p)
{
	for(int i = 0; i < 3; i++)
		if( isnan(p1[i]) || isinf(p1[i]) )
			throw GCMException( GCMException::MESH_EXCEPTION, "NAN or INF in p1");
	for(int i = 0; i < 3; i++)
		if( isnan(p2[i]) || isinf(p2[i]) )
			throw GCMException( GCMException::MESH_EXCEPTION, "NAN or INF in p2");
	for(int i = 0; i < 3; i++)
		if( isnan(p3[i]) || isinf(p3[i]) )
			throw GCMException( GCMException::MESH_EXCEPTION, "NAN or INF in p3");
	for(int i = 0; i < 3; i++)
		if( isnan(p0[i]) || isinf(p0[i]) )
			throw GCMException( GCMException::MESH_EXCEPTION, "NAN or INF in p0");
	for(int i = 0; i < 3; i++)
		if( isnan(v[i]) || isinf(v[i]) )
			throw GCMException( GCMException::MESH_EXCEPTION, "NAN or INF in v");


	bool result;
	// p is point of intersection

	// face normal
	float n[3];
	// face plane parameter
	float d;
	// distance
	float t;

	// Get face normal
	find_border_elem_normal(p1, p2, p3, &n[0], &n[1], &n[2]);

	for(int i = 0; i < 3; i++)
		if( isnan(n[i]) || isinf(n[i]) )
			throw GCMException( GCMException::MESH_EXCEPTION, "NAN or INF in n");

	// If vector is parallel to face - no intersection
	float vn = qm_engine.scalar_product(n[0], n[1], n[2], v[0], v[1], v[2]);
	if( vn == 0 )
		return false;

	// find plane parameter
	d = - qm_engine.scalar_product(n[0], n[1], n[2], p1[0], p1[1], p1[2]);

	// find distance to the plane
	t = - (qm_engine.scalar_product(n[0], n[1], n[2], p0[0], p0[1], p0[2]) + d) / vn;

	// If distance is too big - no intersection
	// If we need opposite direction - no intersection as well
	if( (t < 0) || (t > l) )
		result = false;

	// find point of intersection with the plane
	for(int i = 0; i < 3; i++)
		p[i] = p0[i] + t * v[i];

	// all tests passed - it really intersects
	result = true;

	// check that point is inside triangle
	if( ! qm_engine.same_orientation(p1, p2, p3, p) )
		result = false;
	if( ! qm_engine.same_orientation(p1, p3, p2, p) )
		result = false;
	if( ! qm_engine.same_orientation(p2, p3, p1, p) )
		result = false;

	bool res1;
	float area = fabs( qm_engine.tri_area(p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2], p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]) );
	float areas[3];
	areas[0] = fabs( qm_engine.tri_area(p3[0] - p[0], p3[1] - p[1], p3[2] - p[2], p2[0] - p[0], p2[1] - p[1], p2[2] - p[2]) );
	areas[1] = fabs( qm_engine.tri_area(p1[0] - p[0], p1[1] - p[1], p1[2] - p[2], p2[0] - p[0], p2[1] - p[1], p2[2] - p[2]) );
	areas[2] = fabs( qm_engine.tri_area(p3[0] - p[0], p3[1] - p[1], p3[2] - p[2], p1[0] - p[0], p1[1] - p[1], p1[2] - p[2]) );
	if( fabs(areas[0] + areas[1] + areas[2] - area) > area * 0.0001)
		res1 = false;
//	else if ( (areas[0] == 0) || (areas[1] == 0) || (areas[2] == 0) )
//		res1 = false;
	else
		res1 = true;

//		*logger << "P: " << p[0] << " " << p[1] << " " < p[2];
/*	if(res1 != result) {
		*logger << "Orient res: " << result << " Area res: " << res1 << " with " < (areas[0] + areas[1] + areas[2] > area * 1.00001);
		*logger << "Orients: " << qm_engine.same_orientation(p1, p2, p3, p) << " " 
				<< qm_engine.same_orientation(p1, p3, p2, p) << " "
				< qm_engine.same_orientation(p2, p3, p1, p);
		*logger << "Area: " << area << " Areas: " << areas[0] << " " << areas[1] << " " << areas[2] << " Sum: " < (areas[0] + areas[1] + areas[2]);
		*logger << "P0: " << p1[0] << " " << p1[1] << " " < p1[2];
		*logger << "P1: " << p1[0] << " " << p1[1] << " " < p1[2];
		*logger << "P2: " << p2[0] << " " << p2[1] << " " < p2[2];
		*logger << "P3: " << p3[0] << " " << p3[1] << " " < p3[2];
		*logger << "P: " << p[0] << " " << p[1] << " " < p[2];
//		throw GCMException( GCMException::MESH_EXCEPTION, "Different results for vector_intersects_triangle");
	}*/

	return res1;
};

bool TetrMesh_1stOrder::interpolate_triangle(float *p1, float *p2, float *p3, float *p, float *v1, float *v2, float *v3, float *v, int n)
{
	for(int i = 0; i < 3; i++)
		if( isnan(p[i]) || isinf(p[i]) )
			throw GCMException( GCMException::MESH_EXCEPTION, "NAN or INF in virt node coords");

	float areas[3];
	areas[0] = fabs( qm_engine.tri_area(p3[0] - p[0], p3[1] - p[1], p3[2] - p[2], p2[0] - p[0], p2[1] - p[1], p2[2] - p[2]) );
	areas[1] = fabs( qm_engine.tri_area(p1[0] - p[0], p1[1] - p[1], p1[2] - p[2], p2[0] - p[0], p2[1] - p[1], p2[2] - p[2]) );
	areas[2] = fabs( qm_engine.tri_area(p3[0] - p[0], p3[1] - p[1], p3[2] - p[2], p1[0] - p[0], p1[1] - p[1], p1[2] - p[2]) );

	float l = areas[0] + areas[1] + areas[2];
	float _l1 = areas[0] / l;
	float _l2 = areas[2] / l;
	float _l3 = areas[1] / l;

	// interpolate
	for(int i = 0; i < n; i++)
		v[i] = _l1*v1[i] + _l2*v2[i] + _l3*v3[i];
	return true;

/////////////////////////////

	float n1, n2, n3;
	float p1l[3];
	float p2l[3];
	float p3l[3];
	float pl[3];

	// get face normal
	find_border_elem_normal(p1, p2, p3, &n1, &n2, &n3);

	// FIXME 
	// should this define be moved outside of this function?
	#define sqr(x) (x*x)

	// the new basis:
	// [n1, n2, n3]
	// [n3, 0, -n1]
	// [n1*n2, -(n1^2+n3^2), n2*n3]

	if (n3 != 0.0 || n1 != 0.0)
	{
		gsl_matrix_set(T, 0, 0, n1);
		gsl_matrix_set(T, 0, 1, n2);
		gsl_matrix_set(T, 0, 2, n3);

		pl[0] = sqrt(sqr(n1)+sqr(n3));
		gsl_matrix_set(T, 1, 0, n3/pl[0]);
		gsl_matrix_set(T, 1, 1, 0);
		gsl_matrix_set(T, 1, 2, -n1/pl[0]);

		pl[0] = -(sqr(n1)+sqr(n3));
		pl[1] = sqrt(sqr(n1*n2)+sqr(sqr(n1)+sqr(n3))+sqr(n2*n3));
		gsl_matrix_set(T, 2, 0, n1*n2/pl[1]);
		gsl_matrix_set(T, 2, 1, pl[0]/pl[1]);
		gsl_matrix_set(T, 2, 2, n2*n3/pl[1]);

		// transpose
		gsl_matrix_transpose(T);

		// invert matrix
		int s;
		gsl_linalg_LU_decomp(T, P, &s);
		gsl_linalg_LU_invert(T, P, S);

		// get coordinates in the new basis
		for  (int i = 0; i < 3; i++)
		{

			p1l[i] = 0.0;
			p2l[i] = 0.0;
			p3l[i] = 0.0;
			pl[i] = 0.0;
			for (int j = 0; j < 3; j++)
			{
				p1l[i] += gsl_matrix_get(S, i, j)*p1[j];
				p2l[i] += gsl_matrix_get(S, i, j)*p2[j];
				p3l[i] += gsl_matrix_get(S, i, j)*p3[j];
				pl[i]  += gsl_matrix_get(S, i, j)*p[j];
			}
		}
	}
	else
	{
		p1l[1] = p1[0];
		p1l[2] = p1[2];
		p2l[1] = p2[0];
		p2l[2] = p2[2];
		p3l[1] = p3[0];
		p3l[2] = p3[2];
		pl[1] = p[0];
		pl[2] = p[2];
	}

	// in the new basis we ignore coord #0, because it's the same for all points
	// in face
	
	// get barycentric coordinates for point to be interpolated
	float l1 = ((p2l[2]-p3l[2])*(pl[1]-p3l[1])+(p3l[1]-p2l[1])*(pl[2]-p3l[2]))/((p2l[2]-p3l[2])*(p1l[1]-p3l[1])+(p3l[1]-p2l[1])*(p1l[2]-p3l[2]));
	float l2 = ((p3l[2]-p1l[2])*(pl[1]-p3l[1])+(p1l[1]-p3l[1])*(pl[2]-p3l[2]))/((p2l[2]-p3l[2])*(p1l[1]-p3l[1])+(p3l[1]-p2l[1])*(p1l[2]-p3l[2]));
	float l3 = 1-l2-l1;

	// interpolate
	for(int i = 0; i < n; i++)
		v[i] = l1*v1[i] + l2*v2[i] + l3*v3[i];

	// check if point is inside of face
	return (l1 >= 0.0 && l1 <= 1.0 && l2 >= 0.0 && l2 <= 1.0 && l3 >= 0.0 && l3 <= 1.0); 
};

void TetrMesh_1stOrder::update_current_time(float time_step)
{
	current_time += time_step;
};

void TetrMesh_1stOrder::load_geometry_from_file(string file_name, map<string,string> params){
	enum TYPES {CAS, ELE, GMV, MSH, VTU, ANI3D, UNKNOWN};
	string types_str[] = {"cas", "ele", "gmv", "msh", "vtu", "ani3d"};
	string type = params.count("type") ? params["type"] : "msh";
	int idx = UNKNOWN;
	for (int i = 0; i < UNKNOWN; i++)
		if (type == types_str[i])
		{
			idx = i;
			break;
		}
		switch (idx) {
			// TODO add proper loaders calls
			case CAS:
			{
				load_cas_file(file_name, params.count("id") ? strtol(params["id"].c_str(), NULL, 16) : -1); 
				break;
			}
//			case ELE: ;
//			case GMV: ;
			case MSH:
			{
				load_msh_file(const_cast<char*>(file_name.c_str()));
				break;
			}
			case VTU:
			{
				load_vtu_file(file_name);
				break;
			}
			case ANI3D:
                        {
                                load_ani3d_out_file(const_cast<char*>(file_name.c_str()));
                                break;
                        }
			default: throw GCMException(GCMException::UNIMPLEMENTED_EXCEPTION, "Type "+type+" is not supported");                   
	}
}

void TetrMesh_1stOrder::load_cas_file(string file_name, int zone_id) {
	// TODO get rid of this awesome macro
	#define go_to_end_of_section(f) \
	do { \
		char c = 0; \
		int flag = 1; \
		while (f.good()) \
		{ \
			f >> std::dec >> c; \
			if (c == '(') \
				flag++; \
			else \
				if (c == ')') \
					if (!(--flag)) \
						break; \
		} \
		assert(c == ')'); \
	} while(0) 

	#define check_input(f) \
	do { \
		if (f.eof() || f.fail() || f.bad()) \
			throw GCMException(GCMException::INPUT_EXCEPTION, "Input file is corrupted.");\
	} while (0)

	fstream f;
     
	*logger < "Reading geometry from case file";
	f.open(file_name.c_str());
        
	if (zone_id < 0)
		throw GCMException(GCMException::CONFIG_EXCEPTION, "Invalid zone id");
 
	char c;
	int section, zoneid, fidx, lidx, type, dim;
	int tetrsf = -1, tetrsl = -1;
	int w[5];
	int x;
	vector< set<int> > tetr_vertices;
	float *tmp_nodes;
	int *new_nums;
	
	while (f.good())
	{
		f >> std::dec >> c;
		if (f.eof())
			break;
		if (c == '(')
		{
			f >> std::dec >> section;
			switch (section)
			{
				// ND
				case 2: 
					f >> std::dec >> dim;
					check_input(f);
					assert(dim == 3);
					go_to_end_of_section(f);
					break;
				// NODE
				case  10:
					f >> std::dec >> c;
					check_input(f);
					assert(c == '(');
					f >>  std::hex >> zoneid >> fidx >> lidx >> std::dec >> type >> dim;
					check_input(f);
					// FIXME according to documentation type is 1
					//assert(type == 1); 
					assert(dim == 3);
					go_to_end_of_section(f);
					if (!zoneid)
					{
						tmp_nodes = new float[3*(lidx-fidx+1)];
						new_nums = new int[lidx-fidx+1];
						for (int q = fidx; q <= lidx; q++)
							new_nums[q-fidx] = -1;
						*logger << "Allocted memory for " << (lidx-fidx+1) < " nodes";
					}
					else
					{
						f >> std::dec >> c;
						check_input(f);
						assert(c == '(');
						for (int k = fidx; k <= lidx; k++)
						{
							f >> tmp_nodes[3*(k-fidx)] >> tmp_nodes[3*(k-fidx)+1] >> tmp_nodes[3*(k-fidx)+2];
//							*logger << "NODE " << tmp_nodes[3*(k-fidx)] << " " << " " << tmp_nodes[3*(k-fidx)+1] << " " < tmp_nodes[3*(k-fidx)+2];
							check_input(f);
						}
						go_to_end_of_section(f);						
					}
					go_to_end_of_section(f);
					break;
				// CELL
				case 12:
				{
					f >> std::dec >> c;
					check_input(f);
					assert(c == '(');
					// FIXME we ignore type and element-type since it's not clear
					// if this params are mandatory					
					f >> std::hex >> zoneid >> fidx >> lidx;					
					check_input(f);
					if (zone_id == zoneid)
					{
						tetrsf = fidx;
						tetrsl = lidx;
						tetr_vertices.reserve(lidx-fidx+1);
						for (int k = fidx; k <= lidx; k++)
						{
							tetr_vertices.push_back(set<int>());
						}
						*logger << "Created " << tetr_vertices.size() << " empty tetrs for zone " < zoneid;
					}
					go_to_end_of_section(f);
					// FIXME we do not support cellls' body yet
					go_to_end_of_section(f);
					break;
				}
				// FACES
				case 13:
					f >> std::dec >> c;
					check_input(f);
					assert(c == '(');
					// FIXME we ignore type and element-type since it's not clear
					// if this params are mandatory						
					f >> std::hex >> zoneid >> fidx >> lidx;					
					check_input(f);
					go_to_end_of_section(f);
					if (zoneid)
					{
						f >> std::dec >> c;
						check_input(f);
						assert(c == '(');
						for (int k = fidx; k <= lidx; k++)
						{
							f >> std::hex >> w[0] >> w[1] >> w[2] >> w[3] >> w[4];
							check_input(f);
//							*logger << "FACE " << w[0] << " " << " " << w[1] << " " << w[2] << " " << w[3] << " " < w[4];
							for (int l = 3; l <= 4; l++ )
								if ((w[l] >= tetrsf) && (w[l] <= tetrsl))
									for (int q = 0; q < 3; q++)
									{
										tetr_vertices[w[l]-tetrsf].insert(w[q]);
									}
						}
						go_to_end_of_section(f);
					}
					go_to_end_of_section(f);
					break;					
				default:
					go_to_end_of_section(f);
			}
		}
		else
		{
			throw GCMException(GCMException::INPUT_EXCEPTION, "Unexpected symbol found: "+c);
		}
	}
	
	
	int s = 0;
	for (int i = 0; i < tetrsl-tetrsf+1; i++)
	{
		if (tetr_vertices[i].size() != 4)
			throw GCMException(GCMException::INPUT_EXCEPTION, "Not all vertices read");
		Tetrahedron_1st_order tetr;
		tetr.local_num = i;
		set<int>::iterator iter;
		int q = 0;
		for (iter = tetr_vertices[i].begin(); iter != tetr_vertices[i].end(); iter++)
		{
			int n = (*iter)-1;
//			*logger << "Processing verticle " < n;
			if (new_nums[n] == -1)
			{
				new_nums[n] = s;
				ElasticNode node;
				memcpy(node.coords, tmp_nodes+3*n, 3*sizeof(float));
				node.local_num = s;
				node.elements = NULL;
				node.contact_data = NULL;
				node.local_basis = NULL;
				node.setIsBorder (false);
				node.setContactType (Free);
				node.volume_calculator = NULL;
				node.border_condition = NULL;
				node.contact_condition = NULL;
				node.setPlacement (Local);
				nodes.push_back(node);
				s++;
			}
			tetr.vert[q++] = new_nums[n];
		}
		tetrs.push_back(tetr);
	}
	
	tetr_vertices.clear();
	delete[] tmp_nodes;
	delete[] new_nums;
	f.close();
	*logger << "Loaded " << nodes.size() < " nodes";
	*logger << "Loaded " << tetrs.size() < " tetrs";
}
